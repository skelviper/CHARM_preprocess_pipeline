#!/usr/bin/env python3
"""Run the configured Cutadapt executable only when it is version 4.6."""

import argparse
import os
import subprocess
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a Cutadapt command is required after --")
    return args


def main(argv=None):
    args = parse_args(argv)
    executable = os.path.abspath(args.executable)
    if not os.path.isfile(executable) or not os.access(executable, os.X_OK):
        print("Cutadapt 4.6 executable is unavailable: {}".format(executable), file=sys.stderr)
        return 2
    try:
        version = subprocess.check_output(
            [executable, "--version"], stderr=subprocess.STDOUT, universal_newlines=True
        ).strip()
    except (OSError, subprocess.CalledProcessError) as error:
        print("cannot inspect Cutadapt executable: {}".format(error), file=sys.stderr)
        return 2
    if version != "4.6":
        print(
            "R2 poly(T) trimming requires Cutadapt 4.6, observed {!r}".format(version),
            file=sys.stderr,
        )
        return 2
    os.execv(executable, [executable] + args.command)


if __name__ == "__main__":
    sys.exit(main())
