#!/usr/bin/env python3
"""Discover, freeze, and validate the CHARM Rawdata input contract."""

import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import stat
import sys
import tempfile
import time


SCHEMA_VERSION = "charm_input_contract_v1"
SAMPLE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAFE_CODE_PATTERN = re.compile(r"^C[0-9a-f]{16}$")
FASTQ_SUFFIXES = (".fq.gz", ".fastq.gz")


class InputContractError(ValueError):
    pass


def safe_code_for_sample(sample_name):
    try:
        encoded = sample_name.encode("ascii")
    except UnicodeEncodeError:
        raise InputContractError(
            "sample name must contain only ASCII characters: {!r}".format(sample_name)
        )
    if not SAMPLE_NAME_PATTERN.match(sample_name):
        raise InputContractError(
            "invalid sample name {!r}; expected [A-Za-z0-9][A-Za-z0-9._-]*".format(
                sample_name
            )
        )
    return "C" + hashlib.sha256(encoded).hexdigest()[:16]


def _canonical_json_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def _pretty_json_bytes(value):
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
    ).encode("ascii")


def _payload_sha256(value):
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _strip_line_ending(line):
    return line.rstrip(b"\r\n")


def _first_fastq_record_sha256(path):
    try:
        with gzip.open(path, "rb") as handle:
            record = [handle.readline() for _ in range(4)]
    except (EOFError, OSError) as error:
        raise InputContractError(
            "cannot read gzip FASTQ {!r}: {}".format(path, error)
        )

    if not record[0]:
        raise InputContractError("logically empty FASTQ: {!r}".format(path))
    if any(not line for line in record[1:]):
        raise InputContractError("truncated first FASTQ record: {!r}".format(path))
    if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
        raise InputContractError("invalid first FASTQ record: {!r}".format(path))

    sequence = _strip_line_ending(record[1])
    quality = _strip_line_ending(record[3])
    if not sequence:
        raise InputContractError("empty sequence in first FASTQ record: {!r}".format(path))
    if len(sequence) != len(quality):
        raise InputContractError(
            "sequence/quality length mismatch in first FASTQ record: {!r}".format(path)
        )
    return hashlib.sha256(b"".join(record)).hexdigest()


def _stat_payload(file_stat):
    return {
        "device": int(file_stat.st_dev),
        "inode": int(file_stat.st_ino),
        "mode": int(file_stat.st_mode),
        "mtime_ns": int(file_stat.st_mtime_ns),
        "size": int(file_stat.st_size),
    }


def _logical_path(work_dir, path):
    return os.path.relpath(path, work_dir).replace(os.sep, "/")


def _inspect_fastq(work_dir, path):
    absolute_path = os.path.abspath(path)
    try:
        logical_stat = os.lstat(absolute_path)
    except OSError as error:
        raise InputContractError(
            "cannot inspect FASTQ path {!r}: {}".format(absolute_path, error)
        )

    symlink_target = None
    if stat.S_ISLNK(logical_stat.st_mode):
        try:
            symlink_target = os.readlink(absolute_path)
        except OSError as error:
            raise InputContractError(
                "cannot read FASTQ symlink {!r}: {}".format(absolute_path, error)
            )

    try:
        target_stat = os.stat(absolute_path)
    except OSError as error:
        kind = "broken or looping symlink" if symlink_target is not None else "unreadable path"
        raise InputContractError(
            "{} {!r}: {}".format(kind, absolute_path, error)
        )
    if not stat.S_ISREG(target_stat.st_mode):
        raise InputContractError("FASTQ is not a regular file: {!r}".format(absolute_path))
    if target_stat.st_size == 0:
        raise InputContractError("zero-byte FASTQ: {!r}".format(absolute_path))

    resolved_path = os.path.realpath(absolute_path)
    return {
        "logical_path": _logical_path(work_dir, absolute_path),
        "absolute_path": absolute_path,
        "resolved_path": resolved_path,
        "symlink_target": symlink_target,
        "logical_stat": _stat_payload(logical_stat),
        "target_stat": _stat_payload(target_stat),
        "first_record_sha256": _first_fastq_record_sha256(absolute_path),
    }


def _sample_pair_paths(sample_dir, sample_name):
    return {
        "underscore_numeric": (
            os.path.join(sample_dir, sample_name + "_1.fq.gz"),
            os.path.join(sample_dir, sample_name + "_2.fq.gz"),
        ),
        "illumina_r": (
            os.path.join(sample_dir, sample_name + "_R1.fq.gz"),
            os.path.join(sample_dir, sample_name + "_R2.fq.gz"),
        ),
    }


def _select_pair(sample_dir, sample_name):
    schemes = _sample_pair_paths(sample_dir, sample_name)
    expected_names = {
        os.path.basename(path) for pair in schemes.values() for path in pair
    }
    try:
        names = os.listdir(sample_dir)
    except OSError as error:
        raise InputContractError(
            "cannot enumerate sample directory {!r}: {}".format(sample_dir, error)
        )
    unexpected_fastqs = sorted(
        name
        for name in names
        if not name.startswith(".")
        and name.endswith(FASTQ_SUFFIXES)
        and name not in expected_names
    )
    if unexpected_fastqs:
        raise InputContractError(
            "sample/prefix mismatch or unsupported FASTQ name in {!r}: {}".format(
                sample_dir, ", ".join(unexpected_fastqs)
            )
        )

    present = {
        scheme: tuple(os.path.lexists(path) for path in pair)
        for scheme, pair in schemes.items()
    }
    complete = [scheme for scheme, flags in present.items() if all(flags)]
    any_present = [
        os.path.basename(path)
        for scheme, pair in schemes.items()
        for path, exists in zip(pair, present[scheme])
        if exists
    ]

    if len(complete) == 2:
        raise InputContractError(
            "both supported FASTQ naming schemes are complete for sample {!r}".format(
                sample_name
            )
        )
    selected = complete[0] if len(complete) == 1 else None
    other_scheme_has_files = any(
        any(present[other]) for other in schemes if other != selected
    )
    if selected is None or other_scheme_has_files:
        detail = ", ".join(sorted(any_present)) if any_present else "none"
        raise InputContractError(
            "sample {!r} has a mixed, orphaned, or missing FASTQ pair; present: {}".format(
                sample_name, detail
            )
        )
    return selected, schemes[selected]


def _sample_receipt_payload(contract, sample):
    payload = {
        "schema_version": SCHEMA_VERSION,
        "work_dir": contract["work_dir"],
        "raw_dir": contract["raw_dir"],
        "sample": sample,
    }
    payload["sample_contract_sha256"] = _payload_sha256(payload)
    return payload


def discover_contract(work_dir):
    work_dir = os.path.abspath(work_dir)
    raw_dir = os.path.join(work_dir, "Rawdata")
    if not os.path.isdir(raw_dir):
        raise InputContractError(
            "Rawdata directory is missing or not a directory: {!r}".format(raw_dir)
        )

    try:
        visible_entries = sorted(
            name for name in os.listdir(raw_dir) if not name.startswith(".")
        )
    except OSError as error:
        raise InputContractError("cannot enumerate Rawdata: {}".format(error))
    if not visible_entries:
        raise InputContractError("Rawdata contains no non-hidden sample directories")

    samples = []
    seen_file_keys = {}
    seen_resolved_paths = {}
    seen_safe_codes = {}
    for sample_name in visible_entries:
        safe_code = safe_code_for_sample(sample_name)
        previous_sample = seen_safe_codes.get(safe_code)
        if previous_sample is not None:
            raise InputContractError(
                "safe-code collision for samples {!r} and {!r}: {}".format(
                    previous_sample, sample_name, safe_code
                )
            )
        seen_safe_codes[safe_code] = sample_name

        sample_dir = os.path.join(raw_dir, sample_name)
        if not os.path.isdir(sample_dir):
            raise InputContractError(
                "non-hidden Rawdata entry is not a sample directory: {!r}".format(
                    sample_dir
                )
            )
        naming_scheme, pair = _select_pair(sample_dir, sample_name)
        reads = {}
        for mate, path in zip(("r1", "r2"), pair):
            read = _inspect_fastq(work_dir, path)
            target_stat = read["target_stat"]
            file_key = (target_stat["device"], target_stat["inode"])
            label = "{}:{}".format(sample_name, mate)
            if file_key in seen_file_keys:
                raise InputContractError(
                    "FASTQ target/inode is reused by {} and {}".format(
                        seen_file_keys[file_key], label
                    )
                )
            if read["resolved_path"] in seen_resolved_paths:
                raise InputContractError(
                    "FASTQ resolved target is reused by {} and {}".format(
                        seen_resolved_paths[read["resolved_path"]], label
                    )
                )
            seen_file_keys[file_key] = label
            seen_resolved_paths[read["resolved_path"]] = label
            reads[mate] = read

        receipt_path = "qc/input_contract/samples/{}.json".format(safe_code)
        samples.append(
            {
                "sample_name": sample_name,
                "safe_code": safe_code,
                "naming_scheme": naming_scheme,
                "receipt_path": receipt_path,
                "reads": reads,
            }
        )

    contract = {
        "schema_version": SCHEMA_VERSION,
        "work_dir": work_dir,
        "raw_dir": raw_dir,
        "sample_count": len(samples),
        "samples": samples,
    }
    contract["contract_sha256"] = _payload_sha256(contract)
    return contract


def _contract_core(contract):
    core = dict(contract)
    core.pop("contract_sha256", None)
    return core


def load_contract_file(path):
    try:
        with open(path, "r") as handle:
            contract = json.load(handle)
    except (OSError, ValueError) as error:
        raise InputContractError(
            "cannot load input contract {!r}: {}".format(path, error)
        )
    if contract.get("schema_version") != SCHEMA_VERSION:
        raise InputContractError(
            "unsupported input contract schema in {!r}: {!r}".format(
                path, contract.get("schema_version")
            )
        )
    observed_digest = contract.get("contract_sha256")
    expected_digest = _payload_sha256(_contract_core(contract))
    if observed_digest != expected_digest:
        raise InputContractError(
            "input contract checksum mismatch in {!r}".format(path)
        )

    samples = contract.get("samples")
    if not isinstance(samples, list) or not samples:
        raise InputContractError("input contract has no samples: {!r}".format(path))
    names = [sample.get("sample_name") for sample in samples]
    codes = [sample.get("safe_code") for sample in samples]
    if names != sorted(names) or len(names) != len(set(names)):
        raise InputContractError("input contract sample order/names are invalid")
    if len(codes) != len(set(codes)) or any(
        not isinstance(code, str) or not SAFE_CODE_PATTERN.match(code) for code in codes
    ):
        raise InputContractError("input contract safe codes are invalid or duplicated")
    return contract


def _discovered_cells_tsv(contract):
    columns = [
        "sample_name",
        "safe_code",
        "naming_scheme",
        "r1_logical_path",
        "r2_logical_path",
        "r1_resolved_path",
        "r2_resolved_path",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for sample in contract["samples"]:
        writer.writerow(
            {
                "sample_name": sample["sample_name"],
                "safe_code": sample["safe_code"],
                "naming_scheme": sample["naming_scheme"],
                "r1_logical_path": sample["reads"]["r1"]["logical_path"],
                "r2_logical_path": sample["reads"]["r2"]["logical_path"],
                "r1_resolved_path": sample["reads"]["r1"]["resolved_path"],
                "r2_resolved_path": sample["reads"]["r2"]["resolved_path"],
            }
        )
    return buffer.getvalue().encode("utf-8")


def write_if_changed(path, content):
    try:
        with open(path, "rb") as handle:
            if handle.read() == content:
                return False
    except FileNotFoundError:
        pass
    except OSError as error:
        raise InputContractError("cannot read generated receipt {!r}: {}".format(path, error))

    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".receipt.", dir=parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    return True


def _expected_artifacts(contract):
    artifacts = {
        "current.json": _pretty_json_bytes(contract),
        "discovered_cells.tsv": _discovered_cells_tsv(contract),
    }
    for sample in contract["samples"]:
        relative_path = sample["receipt_path"].split("qc/input_contract/", 1)[1]
        artifacts[relative_path] = _pretty_json_bytes(
            _sample_receipt_payload(contract, sample)
        )
    return artifacts


def _contract_paths(work_dir):
    work_dir = os.path.abspath(work_dir)
    state_root_dir = os.path.join(work_dir, "qc")
    return {
        "work_dir": work_dir,
        "state_root_dir": state_root_dir,
        "state_dir": os.path.join(state_root_dir, "input_contract"),
        "archive_dir": os.path.join(state_root_dir, "input_contract_archive"),
        "lock_path": os.path.join(state_root_dir, ".input_contract.lock"),
    }


def _write_contract_stage(paths, contract):
    stage_dir = tempfile.mkdtemp(
        prefix=".input_contract.stage.", dir=paths["state_root_dir"]
    )
    try:
        for relative_path, content in _expected_artifacts(contract).items():
            artifact_path = os.path.join(stage_dir, relative_path)
            os.makedirs(os.path.dirname(artifact_path), exist_ok=True)
            with open(artifact_path, "wb") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
        return stage_dir
    except BaseException:
        shutil.rmtree(stage_dir, ignore_errors=True)
        raise


def _fsync_directory(path):
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _create_contract_locked(paths):
    state_dir = paths["state_dir"]
    if os.path.lexists(state_dir):
        raise InputContractError(
            "input contract state already exists at {!r}; validate it instead of "
            "recreating it, or use --accept-input-change for an intentional "
            "replacement".format(state_dir)
        )

    contract = discover_contract(paths["work_dir"])
    stage_dir = _write_contract_stage(paths, contract)
    try:
        os.replace(stage_dir, state_dir)
        _fsync_directory(paths["state_root_dir"])
    except BaseException:
        shutil.rmtree(stage_dir, ignore_errors=True)
        raise
    return contract


def create_contract(work_dir):
    """Create the first frozen contract; never replace an existing state."""
    paths = _contract_paths(work_dir)
    os.makedirs(paths["state_root_dir"], exist_ok=True)
    with open(paths["lock_path"], "a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        return _create_contract_locked(paths)


def _validate_frozen_artifacts(work_dir):
    work_dir = os.path.abspath(work_dir)
    contract_path = os.path.join(work_dir, "qc", "input_contract", "current.json")
    contract = load_contract_file(contract_path)
    if contract.get("work_dir") != work_dir:
        raise InputContractError(
            "input contract belongs to {!r}, not {!r}; restore the correct frozen "
            "contract or use --accept-input-change intentionally".format(
                contract.get("work_dir"), work_dir
            )
        )

    state_dir = os.path.dirname(contract_path)
    for relative_path, expected_content in _expected_artifacts(contract).items():
        artifact_path = os.path.join(state_dir, relative_path)
        try:
            with open(artifact_path, "rb") as handle:
                observed_content = handle.read()
        except OSError as error:
            raise InputContractError(
                "missing/unreadable generated input receipt {!r}: {}; restore it "
                "from qc/input_contract_archive or investigate the damaged "
                "run state".format(artifact_path, error)
            )
        if observed_content != expected_content:
            raise InputContractError(
                "generated input receipt drifted: {!r}; restore it from "
                "qc/input_contract_archive or investigate the damaged run "
                "state".format(artifact_path)
            )
    return contract


def _sample_index(contract):
    return {sample["sample_name"]: sample for sample in contract["samples"]}


def _format_change(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=True)


def contract_diff(stored, live):
    """Return deterministic, human-readable differences between two contracts."""
    differences = []
    stored_samples = _sample_index(stored)
    live_samples = _sample_index(live)
    stored_names = set(stored_samples)
    live_names = set(live_samples)

    added = sorted(live_names - stored_names)
    removed = sorted(stored_names - live_names)
    if added:
        differences.append("samples added: {}".format(", ".join(added)))
    if removed:
        differences.append("samples removed: {}".format(", ".join(removed)))

    sample_fields = ("safe_code", "naming_scheme", "receipt_path")
    read_fields = (
        "logical_path",
        "absolute_path",
        "resolved_path",
        "symlink_target",
        "first_record_sha256",
    )
    stat_groups = ("logical_stat", "target_stat")
    stat_fields = ("device", "inode", "mode", "mtime_ns", "size")
    for sample_name in sorted(stored_names & live_names):
        before = stored_samples[sample_name]
        after = live_samples[sample_name]
        for field in sample_fields:
            if before.get(field) != after.get(field):
                differences.append(
                    "sample {} {} changed: {} -> {}".format(
                        sample_name,
                        field,
                        _format_change(before.get(field)),
                        _format_change(after.get(field)),
                    )
                )
        for mate in ("r1", "r2"):
            before_read = before["reads"][mate]
            after_read = after["reads"][mate]
            for field in read_fields:
                if before_read.get(field) != after_read.get(field):
                    differences.append(
                        "sample {} {} {} changed: {} -> {}".format(
                            sample_name,
                            mate,
                            field,
                            _format_change(before_read.get(field)),
                            _format_change(after_read.get(field)),
                        )
                    )
            for group in stat_groups:
                before_stat = before_read[group]
                after_stat = after_read[group]
                for field in stat_fields:
                    if before_stat.get(field) != after_stat.get(field):
                        differences.append(
                            "sample {} {} {}.{} changed: {} -> {}".format(
                                sample_name,
                                mate,
                                group,
                                field,
                                _format_change(before_stat.get(field)),
                                _format_change(after_stat.get(field)),
                            )
                        )

    if not differences and stored.get("contract_sha256") != live.get("contract_sha256"):
        differences.append(
            "contract metadata changed: {} -> {}".format(
                stored.get("contract_sha256"), live.get("contract_sha256")
            )
        )
    return differences


def load_and_validate_contract(work_dir):
    work_dir = os.path.abspath(work_dir)
    contract = _validate_frozen_artifacts(work_dir)
    live_contract = discover_contract(work_dir)
    if live_contract != contract:
        differences = contract_diff(contract, live_contract)
        raise InputContractError(
            "Rawdata differs from the frozen input contract (stored {}, live {}); "
            "the frozen receipts were not changed. Review the differences below, "
            "then run runCHARM.sh --accept-input-change only if the change is "
            "intentional:\n  - {}".format(
                contract.get("contract_sha256"),
                live_contract.get("contract_sha256"),
                "\n  - ".join(differences),
            )
        )
    return contract


def _archive_basename(contract):
    timestamp = time.strftime("%Y%m%dT%H%M%S", time.gmtime())
    nanoseconds = time.time_ns() % 1000000000
    return "{}.{:09d}Z-{}".format(
        timestamp, nanoseconds, contract["contract_sha256"][:12]
    )


def replace_contract(work_dir):
    """Explicitly replace a frozen contract while retaining the previous state."""
    paths = _contract_paths(work_dir)
    os.makedirs(paths["state_root_dir"], exist_ok=True)
    with open(paths["lock_path"], "a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if not os.path.lexists(paths["state_dir"]):
            return _create_contract_locked(paths), None

        stored = _validate_frozen_artifacts(paths["work_dir"])
        live = discover_contract(paths["work_dir"])
        if live == stored:
            return stored, None

        stage_dir = _write_contract_stage(paths, live)
        os.makedirs(paths["archive_dir"], exist_ok=True)
        archive_path = os.path.join(paths["archive_dir"], _archive_basename(stored))
        while os.path.lexists(archive_path):
            archive_path = os.path.join(
                paths["archive_dir"], _archive_basename(stored)
            )

        archived = False
        try:
            os.replace(paths["state_dir"], archive_path)
            archived = True
            os.replace(stage_dir, paths["state_dir"])
            stage_dir = None
            _fsync_directory(paths["state_root_dir"])
            _fsync_directory(paths["archive_dir"])
        except BaseException as error:
            if archived and not os.path.lexists(paths["state_dir"]):
                try:
                    os.replace(archive_path, paths["state_dir"])
                    archived = False
                    _fsync_directory(paths["state_root_dir"])
                except OSError as rollback_error:
                    raise InputContractError(
                        "input-contract replacement failed ({}) and rollback failed "
                        "({}); the prior contract remains at {!r}".format(
                            error, rollback_error, archive_path
                        )
                    )
            raise
        finally:
            if stage_dir is not None:
                shutil.rmtree(stage_dir, ignore_errors=True)
        return live, archive_path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("create", "validate", "replace"):
        command_parser = subparsers.add_parser(command)
        command_parser.add_argument("--work-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        archive_path = None
        if args.command == "create":
            contract = create_contract(args.work_dir)
        elif args.command == "validate":
            contract = load_and_validate_contract(args.work_dir)
        else:
            contract, archive_path = replace_contract(args.work_dir)
    except (InputContractError, OSError) as error:
        print("input contract error: {}".format(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            dict(
                {
                    "contract": os.path.join(
                        contract["work_dir"],
                        "qc",
                        "input_contract",
                        "current.json",
                    ),
                    "contract_sha256": contract["contract_sha256"],
                    "sample_count": contract["sample_count"],
                },
                **({"archived_contract": archive_path} if archive_path else {})
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
