#!/usr/bin/env python3
"""Resolve the launcher and Snakemake to one effective CHARM work directory."""

import argparse
import ast
import os
import re
import sys

try:
    import yaml
except ImportError as error:
    yaml = None
    YAML_IMPORT_ERROR = error
else:
    YAML_IMPORT_ERROR = None


CONFIG_KEY_PATTERN = re.compile(r"^[A-Za-z_]\w*$")
RUN_OWNED_OPTIONS = (
    "--directory",
    "--profile",
    "--snakefile",
    "-d",
    "-s",
)
class RunConfigError(ValueError):
    pass


def _load_config(path):
    if yaml is None:
        raise RunConfigError(
            "PyYAML is required to resolve config before launch ({}); run the "
            "launcher in the pinned charm environment".format(YAML_IMPORT_ERROR)
        )
    try:
        with open(path, "r") as handle:
            value = yaml.safe_load(handle)
    except (OSError, ValueError, yaml.YAMLError) as error:
        raise RunConfigError("cannot load config file {!r}: {}".format(path, error))
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise RunConfigError("config file must contain a mapping: {!r}".format(path))
    return value


def _update_config(target, update):
    for key, value in update.items():
        if isinstance(value, dict):
            existing = target.get(key)
            if not isinstance(existing, dict):
                existing = {}
            target[key] = _update_config(existing, value)
        else:
            target[key] = value
    return target


def _parse_config_value(value):
    for parser in (int, float, ast.literal_eval, str):
        try:
            parsed = parser(value)
        except (ValueError, SyntaxError):
            continue
        if not callable(parsed):
            return parsed
    return value


def _option_is(argument, long_name, short_name=None):
    if argument == long_name or argument.startswith(long_name + "="):
        return True
    if short_name and (argument == short_name or argument.startswith(short_name)):
        return True
    return False


def _consume_values(arguments, start, attached, option):
    values = [] if attached is None else [attached]
    index = start
    while index < len(arguments) and not arguments[index].startswith("-"):
        values.append(arguments[index])
        index += 1
    if not values:
        raise RunConfigError("{} requires at least one value".format(option))
    return values, index


def parse_snakemake_arguments(arguments, invocation_dir):
    """Parse config-bearing arguments and normalize config-file paths."""
    normalized = []
    configfiles = None
    config_entries = None
    index = 0
    while index < len(arguments):
        argument = arguments[index]
        if argument == "--":
            normalized.extend(arguments[index:])
            break

        if any(_option_is(argument, option) for option in RUN_OWNED_OPTIONS[:3]) or (
            argument.startswith("-d") and argument != "--dry-run" and argument != "--dryrun"
        ) or (argument.startswith("-s") and not argument.startswith("--")):
            raise RunConfigError(
                "runCHARM.sh owns --snakefile and --directory and does not accept "
                "--profile; configure the pipeline through --configfile or "
                "--config work_dir=..."
            )

        if argument in ("--configfile", "--configfiles") or argument.startswith(
            ("--configfile=", "--configfiles=")
        ):
            if configfiles is not None:
                raise RunConfigError(
                    "repeat --configfile is ambiguous in Snakemake 5.20.1; provide "
                    "all files after one --configfile option"
                )
            attached = argument.split("=", 1)[1] if "=" in argument else None
            values, index = _consume_values(
                arguments, index + 1, attached, "--configfile"
            )
            configfiles = [
                path if os.path.isabs(path) else os.path.join(invocation_dir, path)
                for path in values
            ]
            configfiles = [os.path.abspath(path) for path in configfiles]
            normalized.append("--configfile")
            normalized.extend(configfiles)
            continue

        is_config = argument in ("--config", "-C") or argument.startswith(
            "--config="
        ) or (argument.startswith("-C") and argument != "-C")
        if is_config:
            if config_entries is not None:
                raise RunConfigError(
                    "repeat --config is ambiguous in Snakemake 5.20.1; provide all "
                    "KEY=VALUE entries after one --config option"
                )
            if argument.startswith("--config="):
                attached = argument.split("=", 1)[1]
            elif argument.startswith("-C") and argument != "-C":
                attached = argument[2:]
            else:
                attached = None
            values, index = _consume_values(arguments, index + 1, attached, "--config")
            for entry in values:
                if "=" not in entry:
                    raise RunConfigError(
                        "invalid --config value {!r}; values must be KEY=VALUE and "
                        "targets must precede --config".format(entry)
                    )
                key = entry.split("=", 1)[0]
                if not CONFIG_KEY_PATTERN.match(key):
                    raise RunConfigError(
                        "invalid --config key {!r}; expected a Python identifier".format(
                            key
                        )
                    )
            config_entries = values
            normalized.append("--config")
            normalized.extend(values)
            continue

        normalized.append(argument)
        index += 1

    return normalized, configfiles or [], config_entries or []


def resolve(pipeline_dir, invocation_dir, arguments):
    pipeline_dir = os.path.abspath(pipeline_dir)
    invocation_dir = os.path.abspath(invocation_dir)
    normalized, configfiles, config_entries = parse_snakemake_arguments(
        list(arguments), invocation_dir
    )

    config_path = os.path.join(pipeline_dir, "config.yaml")
    config = _load_config(config_path)
    for path in configfiles:
        _update_config(config, _load_config(path))
    for entry in config_entries:
        key, value = entry.split("=", 1)
        config[key] = _parse_config_value(value)

    work_dir = config.get("work_dir", "..")
    if not isinstance(work_dir, str):
        raise RunConfigError(
            "effective work_dir must be a string, observed {!r}".format(work_dir)
        )
    work_dir = work_dir.replace("\r", "").strip()
    if not work_dir:
        raise RunConfigError("effective work_dir is empty")
    if not os.path.isabs(work_dir):
        work_dir = os.path.join(pipeline_dir, work_dir)
    work_dir = os.path.abspath(os.path.normpath(work_dir))
    if not os.path.isdir(work_dir):
        raise RunConfigError(
            "effective work_dir is missing or not a directory: {!r}".format(work_dir)
        )
    return work_dir, normalized


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline-dir", required=True)
    parser.add_argument("--invocation-dir", required=True)
    parser.add_argument("snakemake_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.snakemake_args[:1] == ["--"]:
        args.snakemake_args = args.snakemake_args[1:]
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        work_dir, normalized = resolve(
            args.pipeline_dir,
            args.invocation_dir,
            args.snakemake_args,
        )
    except RunConfigError as error:
        print("run config error: {}".format(error), file=sys.stderr)
        return 2
    values = [work_dir] + normalized
    sys.stdout.buffer.write(b"\0".join(value.encode("utf-8") for value in values) + b"\0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
