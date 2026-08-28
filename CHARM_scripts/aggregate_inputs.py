#!/usr/bin/env python3
"""Aggregate large Snakemake input lists without expanding them into argv."""

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


COPY_BLOCK_SIZE = 8 * 1024 * 1024


class AggregateInputError(RuntimeError):
    pass


def _validated_inputs(input_paths):
    paths = [Path(str(path)) for path in input_paths]
    if not paths:
        raise AggregateInputError("at least one input file is required")
    for path in paths:
        if not path.is_file():
            raise AggregateInputError("input is not a regular file: {}".format(path))
        if path.stat().st_size == 0:
            raise AggregateInputError("input is empty: {}".format(path))
    return paths


def _temporary_path(output, suffix=""):
    output = Path(str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix="." + output.name + ".", suffix=suffix, dir=str(output.parent)
    )
    os.close(descriptor)
    return Path(name)


def concatenate_gzip_members(input_paths, output_path):
    """Byte-concatenate gzip members in input order and publish atomically."""
    paths = _validated_inputs(input_paths)
    output = Path(str(output_path))
    temporary = _temporary_path(output)
    try:
        with temporary.open("wb") as destination:
            for path in paths:
                with path.open("rb") as source:
                    shutil.copyfileobj(source, destination, COPY_BLOCK_SIZE)
        if temporary.stat().st_size == 0:
            raise AggregateInputError("concatenated gzip output is empty")
        os.replace(str(temporary), str(output))
    finally:
        if temporary.exists():
            temporary.unlink()


def _stderr_text(handle):
    handle.seek(0)
    return handle.read().decode("utf-8", errors="replace").strip()


def _stop_processes(processes):
    for process in processes:
        if process.poll() is None:
            process.terminate()
    for process in processes:
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def combine_fragment_files(input_paths, output_path, index_path, threads):
    """Stream gzip inputs through the production sort/uniq/bgzip/tabix chain."""
    paths = _validated_inputs(input_paths)
    output = Path(str(output_path))
    index = Path(str(index_path))
    output.parent.mkdir(parents=True, exist_ok=True)
    index.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(output, suffix=".bgz")
    temporary_index = Path(str(temporary) + ".tbi")
    processes = []

    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    thread_count = max(1, int(threads))

    try:
        with tempfile.TemporaryFile() as sort_stderr, tempfile.TemporaryFile() as uniq_stderr, tempfile.TemporaryFile() as bgzip_stderr:
            sort_process = subprocess.Popen(
                ["sort", "-k1,1", "-k2,2n"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=sort_stderr,
                env=environment,
            )
            processes.append(sort_process)
            uniq_process = subprocess.Popen(
                ["uniq"],
                stdin=sort_process.stdout,
                stdout=subprocess.PIPE,
                stderr=uniq_stderr,
                env=environment,
            )
            processes.append(uniq_process)
            sort_process.stdout.close()

            with temporary.open("wb") as destination:
                bgzip_process = subprocess.Popen(
                    ["bgzip", "-@", str(thread_count)],
                    stdin=uniq_process.stdout,
                    stdout=destination,
                    stderr=bgzip_stderr,
                    env=environment,
                )
                processes.append(bgzip_process)
                uniq_process.stdout.close()
                try:
                    for path in paths:
                        with gzip.open(str(path), "rb") as source:
                            shutil.copyfileobj(
                                source, sort_process.stdin, COPY_BLOCK_SIZE
                            )
                finally:
                    sort_process.stdin.close()

                return_codes = {
                    "bgzip": bgzip_process.wait(),
                    "uniq": uniq_process.wait(),
                    "sort": sort_process.wait(),
                }
                errors = {
                    "sort": _stderr_text(sort_stderr),
                    "uniq": _stderr_text(uniq_stderr),
                    "bgzip": _stderr_text(bgzip_stderr),
                }
                failures = [
                    "{} exited {}{}".format(
                        name,
                        return_codes[name],
                        ": " + errors[name] if errors[name] else "",
                    )
                    for name in ("sort", "uniq", "bgzip")
                    if return_codes[name] != 0
                ]
                if failures:
                    raise AggregateInputError("; ".join(failures))

        completed = subprocess.run(
            ["tabix", "-f", "-p", "bed", str(temporary)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )
        if completed.returncode != 0:
            message = completed.stderr.decode("utf-8", errors="replace").strip()
            raise AggregateInputError(
                "tabix exited {}{}".format(
                    completed.returncode, ": " + message if message else ""
                )
            )
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise AggregateInputError("combined fragment output is empty")
        if not temporary_index.is_file() or temporary_index.stat().st_size == 0:
            raise AggregateInputError("tabix index is empty")

        os.replace(str(temporary), str(output))
        os.replace(str(temporary_index), str(index))
    except Exception:
        _stop_processes(processes)
        raise
    finally:
        for path in (temporary, temporary_index):
            if path.exists():
                path.unlink()
