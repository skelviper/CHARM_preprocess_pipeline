#!/usr/bin/env python3

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import yaml


PIPELINE = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PIPELINE / "CHARM_scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import resolve_run_config


def write_pair(work_dir, sample):
    sample_dir = Path(work_dir) / "Rawdata" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    for mate in (1, 2):
        with gzip.open(sample_dir / "{}_{}.fq.gz".format(sample, mate), "wt") as handle:
            handle.write("@{}_{}\nACGT\n+\nIIII\n".format(sample, mate))


def state_snapshot(work_dir):
    state_dir = Path(work_dir) / "qc" / "input_contract"
    return {
        str(path.relative_to(state_dir)): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in sorted(state_dir.rglob("*"))
        if path.is_file()
    }


def make_launcher_fixture(tmp_path, work_dir):
    pipeline = tmp_path / "pipeline"
    scripts = pipeline / "CHARM_scripts"
    scripts.mkdir(parents=True)
    shutil.copy2(PIPELINE / "runCHARM.sh", pipeline / "runCHARM.sh")
    shutil.copy2(SCRIPT_DIR / "input_contract.py", scripts / "input_contract.py")
    shutil.copy2(SCRIPT_DIR / "resolve_run_config.py", scripts / "resolve_run_config.py")
    (pipeline / "CHARM.smk").write_text("rule all:\n    shell: 'true'\n")
    (pipeline / "config.yaml").write_text(
        yaml.safe_dump({"work_dir": str(work_dir)}, sort_keys=False)
    )
    fake_snakemake = pipeline / "fake_snakemake.py"
    fake_snakemake.write_text(
        """#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path
capture = os.environ.get("FAKE_SNAKEMAKE_CAPTURE")
if capture:
    Path(capture).write_text(json.dumps({
        "args": sys.argv[1:],
        "cwd": os.getcwd(),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV"),
        "effective_work_dir": os.environ.get("CHARM_EFFECTIVE_WORK_DIR"),
        "tmpdir": os.environ.get("TMPDIR"),
    }, sort_keys=True))
"""
    )
    fake_snakemake.chmod(0o755)
    return pipeline, fake_snakemake


def launcher_environment(fake_snakemake):
    environment = os.environ.copy()
    environment.pop("TMPDIR", None)
    environment.update(
        {
            "CONDA_DEFAULT_ENV": "charm",
            "CHARM_PYTHON_BIN": sys.executable,
            "CHARM_SNAKEMAKE_BIN": str(fake_snakemake),
        }
    )
    return environment


def run_launcher(pipeline, environment, *arguments):
    return subprocess.run(
        [str(pipeline / "runCHARM.sh")] + list(arguments),
        cwd=str(pipeline.parent),
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def test_resolver_merges_configfiles_and_cli_work_dir(tmp_path):
    pipeline = tmp_path / "pipeline"
    pipeline.mkdir()
    default_work = tmp_path / "default"
    file_work = tmp_path / "from-file"
    cli_work = tmp_path / "from-cli"
    for path in (default_work, file_work, cli_work):
        path.mkdir()
    (pipeline / "config.yaml").write_text(
        yaml.safe_dump({"work_dir": str(default_work), "mapping": {"mode": "per_cell"}})
    )
    override = tmp_path / "override.yaml"
    override.write_text(yaml.safe_dump({"work_dir": str(file_work)}))

    work_dir, normalized = resolve_run_config.resolve(
        str(pipeline),
        str(tmp_path),
        ["target", "--configfile", override.name, "--dry-run"],
    )
    assert work_dir == str(file_work)
    assert normalized == ["target", "--configfile", str(override), "--dry-run"]

    work_dir, normalized = resolve_run_config.resolve(
        str(pipeline),
        str(tmp_path),
        ["target", "--config", "work_dir={}".format(cli_work), "--dry-run"],
    )
    assert work_dir == str(cli_work)
    assert normalized[-3:] == [
        "--config",
        "work_dir={}".format(cli_work),
        "--dry-run",
    ]

def test_launcher_freezes_validates_and_explicitly_archives_input_change(tmp_path):
    work_dir = tmp_path / "workspace"
    work_dir.mkdir()
    write_pair(work_dir, "Alpha")
    pipeline, fake_snakemake = make_launcher_fixture(tmp_path, work_dir)
    capture = tmp_path / "snakemake.json"
    environment = launcher_environment(fake_snakemake)
    environment["FAKE_SNAKEMAKE_CAPTURE"] = str(capture)

    first = run_launcher(pipeline, environment, "--dry-run")
    assert first.returncode == 0, first.stderr
    assert (work_dir / "qc" / "logs").is_dir()
    assert not (pipeline / "qc").exists()
    before = state_snapshot(work_dir)
    second = run_launcher(pipeline, environment, "--dry-run")
    assert second.returncode == 0, second.stderr
    assert state_snapshot(work_dir) == before

    write_pair(work_dir, "Beta")
    drift = run_launcher(pipeline, environment, "--dry-run")
    assert drift.returncode == 2
    assert "samples added: Beta" in drift.stderr
    assert state_snapshot(work_dir) == before

    accepted = run_launcher(
        pipeline, environment, "--accept-input-change", "--dry-run"
    )
    assert accepted.returncode == 0, accepted.stderr
    current = json.loads(
        (work_dir / "qc" / "input_contract" / "current.json").read_text()
    )
    assert [sample["sample_name"] for sample in current["samples"]] == ["Alpha", "Beta"]
    archives = list((work_dir / "qc" / "input_contract_archive").iterdir())
    assert len(archives) == 1
    archived = json.loads((archives[0] / "current.json").read_text())
    assert [sample["sample_name"] for sample in archived["samples"]] == ["Alpha"]
    observed_args = json.loads(capture.read_text())["args"]
    assert "--accept-input-change" not in observed_args
    assert observed_args == [
        "--cluster",
        "sbatch --qos=high -w node03 --output=qc/logs/slurm-%j.out "
        "--cpus-per-task={threads} -t 7-00:00 -J CHARM!",
        "--jobs",
        "1024",
        "--resources",
        "star_slots=1",
        "count_slots=1",
        "--rerun-incomplete",
        "-s",
        str(pipeline / "CHARM.smk"),
        "--dry-run",
    ]
    assert "--keep-going" not in observed_args
    assert "--profile" not in observed_args
    observed = json.loads(capture.read_text())
    assert observed["conda_default_env"] == "charm"
    assert observed["tmpdir"] == str(work_dir / "tmp")
    assert not (work_dir / "tmp").exists()


def test_launcher_uses_configfile_workspace_for_contract_and_dag(tmp_path):
    default_work = tmp_path / "default-work"
    override_work = tmp_path / "override-work"
    default_work.mkdir()
    override_work.mkdir()
    write_pair(override_work, "Beta")
    pipeline, fake_snakemake = make_launcher_fixture(tmp_path, default_work)
    override = tmp_path / "override.yaml"
    override.write_text(yaml.safe_dump({"work_dir": str(override_work)}))
    capture = tmp_path / "effective.json"
    environment = launcher_environment(fake_snakemake)
    environment["FAKE_SNAKEMAKE_CAPTURE"] = str(capture)

    result = run_launcher(
        pipeline, environment, "target", "--configfile", str(override), "--dry-run"
    )
    assert result.returncode == 0, result.stderr
    observed = json.loads(capture.read_text())
    assert observed["cwd"] == str(override_work)
    assert observed["effective_work_dir"] == str(override_work)
    assert observed["tmpdir"] == str(override_work / "tmp")
    assert str(override) in observed["args"]
    assert (override_work / "qc" / "input_contract" / "current.json").exists()
    assert not (default_work / "qc").exists()


def test_launcher_preserves_explicit_tmpdir(tmp_path):
    work_dir = tmp_path / "workspace"
    work_dir.mkdir()
    write_pair(work_dir, "Alpha")
    pipeline, fake_snakemake = make_launcher_fixture(tmp_path, work_dir)
    capture = tmp_path / "explicit-tmp.json"
    explicit_tmp = tmp_path / "scratch" / "charm"
    environment = launcher_environment(fake_snakemake)
    environment["FAKE_SNAKEMAKE_CAPTURE"] = str(capture)
    environment["TMPDIR"] = str(explicit_tmp)

    result = run_launcher(pipeline, environment, "--dry-run")
    assert result.returncode == 0, result.stderr
    observed = json.loads(capture.read_text())
    assert observed["tmpdir"] == str(explicit_tmp)
    assert explicit_tmp.is_dir()
    assert not (work_dir / "tmp").exists()


def test_launcher_reenters_existing_charm_environment(tmp_path):
    work_dir = tmp_path / "workspace"
    work_dir.mkdir()
    write_pair(work_dir, "Alpha")
    pipeline, fake_snakemake = make_launcher_fixture(tmp_path, work_dir)
    capture = tmp_path / "reentered.json"
    fake_conda = tmp_path / "fake-conda"
    fake_conda.write_text(
        """#!/usr/bin/env bash
set -euo pipefail
test "$1" = run
shift
test "$1" = --no-capture-output
shift
test "$1" = -n
shift
export CONDA_DEFAULT_ENV="$1"
shift
exec "$@"
"""
    )
    fake_conda.chmod(0o755)
    environment = launcher_environment(fake_snakemake)
    environment.pop("CONDA_DEFAULT_ENV")
    environment["CHARM_CONDA_BIN"] = str(fake_conda)
    environment["FAKE_SNAKEMAKE_CAPTURE"] = str(capture)

    result = run_launcher(pipeline, environment, "--dry-run")
    assert result.returncode == 0, result.stderr
    observed = json.loads(capture.read_text())
    assert observed["conda_default_env"] == "charm"
    assert observed["tmpdir"] == str(work_dir / "tmp")


def test_launcher_rejects_owned_directory_before_qc_mutation(tmp_path):
    work_dir = tmp_path / "workspace"
    work_dir.mkdir()
    write_pair(work_dir, "Alpha")
    pipeline, fake_snakemake = make_launcher_fixture(tmp_path, work_dir)
    environment = launcher_environment(fake_snakemake)

    result = run_launcher(pipeline, environment, "--directory", str(tmp_path))
    assert result.returncode == 2
    assert "runCHARM.sh owns" in result.stderr
    assert not (work_dir / "qc").exists()
