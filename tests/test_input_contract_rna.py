import gzip
import json
import os
import shutil
import sys
from pathlib import Path

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "CHARM_scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import input_contract
from input_contract import InputContractError
from rename_count_matrix import rename_count_matrix
from tag_fastq_cell import tag_fastq


def write_fastq(path, read_id="read", sequence="ACGT", quality="IIII"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("@{}\n{}\n+\n{}\n".format(read_id, sequence, quality))


def add_pair(work_dir, sample, scheme="numeric"):
    sample_dir = Path(work_dir) / "Rawdata" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    if scheme == "numeric":
        names = ("{}_1.fq.gz".format(sample), "{}_2.fq.gz".format(sample))
    else:
        names = ("{}_R1.fq.gz".format(sample), "{}_R2.fq.gz".format(sample))
    write_fastq(sample_dir / names[0], read_id="{}_r1".format(sample))
    write_fastq(sample_dir / names[1], read_id="{}_r2".format(sample))
    return sample_dir, names


def sample_map(contract):
    return {sample["sample_name"]: sample for sample in contract["samples"]}


def state_snapshot(work_dir):
    state_dir = Path(work_dir) / "qc" / "input_contract"
    return {
        str(path.relative_to(state_dir)): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in sorted(state_dir.rglob("*"))
        if path.is_file()
    }


def test_discovers_both_schemes_and_preserves_unchanged_receipt_mtime(tmp_path):
    add_pair(tmp_path, "Beta-2", scheme="illumina")
    add_pair(tmp_path, "Alpha_one", scheme="numeric")

    contract = input_contract.create_contract(str(tmp_path))
    assert [sample["sample_name"] for sample in contract["samples"]] == [
        "Alpha_one",
        "Beta-2",
    ]
    info = sample_map(contract)
    assert info["Alpha_one"]["naming_scheme"] == "underscore_numeric"
    assert info["Beta-2"]["naming_scheme"] == "illumina_r"
    assert info["Alpha_one"]["safe_code"] == input_contract.safe_code_for_sample(
        "Alpha_one"
    )

    state_dir = tmp_path / "qc" / "input_contract"
    artifacts = [
        state_dir / "current.json",
        state_dir / "discovered_cells.tsv",
    ] + [tmp_path / sample["receipt_path"] for sample in contract["samples"]]
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in artifacts}
    second = input_contract.load_and_validate_contract(str(tmp_path))
    after = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in artifacts}
    assert second == contract
    assert after == before
    with pytest.raises(InputContractError, match="already exists"):
        input_contract.create_contract(str(tmp_path))
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in artifacts
    } == before


def test_both_complete_naming_schemes_fail(tmp_path):
    sample_dir, _ = add_pair(tmp_path, "Alpha", scheme="numeric")
    write_fastq(sample_dir / "Alpha_R1.fq.gz")
    write_fastq(sample_dir / "Alpha_R2.fq.gz")
    with pytest.raises(InputContractError, match="both supported"):
        input_contract.discover_contract(str(tmp_path))


@pytest.mark.parametrize(
    "present",
    [
        ("Alpha_1.fq.gz",),
        ("Alpha_1.fq.gz", "Alpha_R2.fq.gz"),
        ("Alpha_R1.fq.gz",),
    ],
)
def test_orphan_and_mixed_pairs_fail(tmp_path, present):
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    for name in present:
        write_fastq(sample_dir / name)
    with pytest.raises(InputContractError, match="mixed, orphaned, or missing"):
        input_contract.discover_contract(str(tmp_path))


def test_zero_byte_fastq_fails(tmp_path):
    sample_dir, names = add_pair(tmp_path, "Alpha")
    (sample_dir / names[0]).write_bytes(b"")
    with pytest.raises(InputContractError, match="zero-byte FASTQ"):
        input_contract.discover_contract(str(tmp_path))


def test_logically_empty_gzip_fastq_fails(tmp_path):
    sample_dir, names = add_pair(tmp_path, "Alpha")
    with gzip.open(sample_dir / names[0], "wb"):
        pass
    with pytest.raises(InputContractError, match="logically empty FASTQ"):
        input_contract.discover_contract(str(tmp_path))


def test_broken_symlink_fails(tmp_path):
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    os.symlink("missing.fq.gz", sample_dir / "Alpha_1.fq.gz")
    write_fastq(sample_dir / "Alpha_2.fq.gz")
    with pytest.raises(InputContractError, match="broken or looping symlink"):
        input_contract.discover_contract(str(tmp_path))


def test_non_regular_fastq_fails(tmp_path):
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    (sample_dir / "Alpha_1.fq.gz").mkdir()
    write_fastq(sample_dir / "Alpha_2.fq.gz")
    with pytest.raises(InputContractError, match="not a regular file"):
        input_contract.discover_contract(str(tmp_path))


def test_prefix_mismatch_fails(tmp_path):
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    write_fastq(sample_dir / "Other_1.fq.gz")
    write_fastq(sample_dir / "Other_2.fq.gz")
    with pytest.raises(InputContractError, match="sample/prefix mismatch"):
        input_contract.discover_contract(str(tmp_path))


@pytest.mark.parametrize("sample", ["has space", "bad$name", "nonascii_\u03b1", "_starts_bad"])
def test_illegal_sample_names_fail(tmp_path, sample):
    add_pair(tmp_path, sample)
    with pytest.raises(InputContractError, match="sample name"):
        input_contract.discover_contract(str(tmp_path))


def test_r1_r2_same_inode_fails(tmp_path):
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    write_fastq(sample_dir / "Alpha_1.fq.gz")
    os.link(sample_dir / "Alpha_1.fq.gz", sample_dir / "Alpha_2.fq.gz")
    with pytest.raises(InputContractError, match="target/inode is reused"):
        input_contract.discover_contract(str(tmp_path))


def test_cross_sample_duplicate_target_fails(tmp_path):
    first_dir, _ = add_pair(tmp_path, "Alpha")
    second_dir = tmp_path / "Rawdata" / "Beta"
    second_dir.mkdir(parents=True)
    os.symlink(
        os.path.relpath(first_dir / "Alpha_1.fq.gz", second_dir),
        second_dir / "Beta_1.fq.gz",
    )
    write_fastq(second_dir / "Beta_2.fq.gz")
    with pytest.raises(InputContractError, match="target/inode is reused"):
        input_contract.discover_contract(str(tmp_path))


def test_discovery_is_stable_under_directory_enumeration_order(tmp_path, monkeypatch):
    add_pair(tmp_path, "Zulu")
    add_pair(tmp_path, "Alpha")
    expected = input_contract.discover_contract(str(tmp_path))
    original_listdir = os.listdir

    def reversed_listdir(path):
        return list(reversed(original_listdir(path)))

    monkeypatch.setattr(input_contract.os, "listdir", reversed_listdir)
    assert input_contract.discover_contract(str(tmp_path)) == expected


def test_adding_sample_fails_without_mutation_then_explicitly_archives_and_replaces(
    tmp_path,
):
    add_pair(tmp_path, "Alpha")
    first = input_contract.create_contract(str(tmp_path))
    alpha = sample_map(first)["Alpha"]
    alpha_receipt = tmp_path / alpha["receipt_path"]
    receipt_mtime = alpha_receipt.stat().st_mtime_ns
    before = state_snapshot(tmp_path)

    add_pair(tmp_path, "Beta")
    with pytest.raises(InputContractError, match="samples added: Beta"):
        input_contract.load_and_validate_contract(str(tmp_path))
    assert state_snapshot(tmp_path) == before

    second, archive_path = input_contract.replace_contract(str(tmp_path))
    assert sample_map(second)["Alpha"]["safe_code"] == alpha["safe_code"]
    assert archive_path is not None
    archive_path = Path(archive_path)
    assert archive_path.parent == tmp_path / "qc" / "input_contract_archive"
    assert (archive_path / "current.json").read_bytes() == before["current.json"][0]
    archived_receipt = archive_path / "samples" / alpha_receipt.name
    assert archived_receipt.stat().st_mtime_ns == receipt_mtime
    assert input_contract.load_and_validate_contract(str(tmp_path)) == second


def test_deleted_sample_fails_and_preserves_frozen_state(tmp_path):
    add_pair(tmp_path, "Alpha")
    add_pair(tmp_path, "Beta")
    input_contract.create_contract(str(tmp_path))
    before = state_snapshot(tmp_path)

    shutil.rmtree(str(tmp_path / "Rawdata" / "Beta"))
    with pytest.raises(InputContractError, match="samples removed: Beta"):
        input_contract.load_and_validate_contract(str(tmp_path))
    assert state_snapshot(tmp_path) == before


def test_symlink_retarget_fails_with_diff_and_preserves_frozen_state(tmp_path):
    delivery = tmp_path / "delivery"
    write_fastq(delivery / "first.R1.fq.gz", read_id="first_r1")
    write_fastq(delivery / "second.R1.fq.gz", read_id="second_r1")
    write_fastq(delivery / "R2.fq.gz", read_id="r2")
    sample_dir = tmp_path / "Rawdata" / "Alpha"
    sample_dir.mkdir(parents=True)
    r1_link = sample_dir / "Alpha_1.fq.gz"
    r2_link = sample_dir / "Alpha_2.fq.gz"
    os.symlink(os.path.relpath(delivery / "first.R1.fq.gz", sample_dir), r1_link)
    os.symlink(os.path.relpath(delivery / "R2.fq.gz", sample_dir), r2_link)
    input_contract.create_contract(str(tmp_path))
    before = state_snapshot(tmp_path)

    r1_link.unlink()
    os.symlink(os.path.relpath(delivery / "second.R1.fq.gz", sample_dir), r1_link)
    with pytest.raises(InputContractError, match="resolved_path changed"):
        input_contract.load_and_validate_contract(str(tmp_path))
    assert state_snapshot(tmp_path) == before


def test_size_and_stat_drift_fail_with_diff_and_preserve_frozen_state(tmp_path):
    sample_dir, names = add_pair(tmp_path, "Alpha")
    r1 = sample_dir / names[0]
    input_contract.create_contract(str(tmp_path))
    before = state_snapshot(tmp_path)

    write_fastq(r1, read_id="changed", sequence="ACGTACGTACGT", quality="I" * 12)
    with pytest.raises(InputContractError, match=r"target_stat\.(size|mtime_ns) changed"):
        input_contract.load_and_validate_contract(str(tmp_path))
    assert state_snapshot(tmp_path) == before

    input_contract.replace_contract(str(tmp_path))
    replaced = state_snapshot(tmp_path)
    observed = r1.stat()
    os.utime(r1, ns=(observed.st_atime_ns, observed.st_mtime_ns + 1000000000))
    with pytest.raises(InputContractError, match=r"target_stat\.mtime_ns changed"):
        input_contract.load_and_validate_contract(str(tmp_path))
    assert state_snapshot(tmp_path) == replaced


def test_failed_explicit_replacement_leaves_previous_contract_intact(tmp_path):
    add_pair(tmp_path, "Alpha")
    input_contract.create_contract(str(tmp_path))
    before = state_snapshot(tmp_path)

    broken_dir = tmp_path / "Rawdata" / "Broken"
    broken_dir.mkdir()
    write_fastq(broken_dir / "Broken_1.fq.gz")
    with pytest.raises(InputContractError, match="mixed, orphaned, or missing"):
        input_contract.replace_contract(str(tmp_path))

    assert state_snapshot(tmp_path) == before
    archive_root = tmp_path / "qc" / "input_contract_archive"
    assert not archive_root.exists() or not list(archive_root.iterdir())


def test_failed_publish_rolls_back_exact_previous_contract(tmp_path, monkeypatch):
    add_pair(tmp_path, "Alpha")
    input_contract.create_contract(str(tmp_path))
    before = state_snapshot(tmp_path)
    add_pair(tmp_path, "Beta")
    real_replace = input_contract.os.replace
    calls = []

    def fail_second_replace(source, destination):
        calls.append((source, destination))
        if len(calls) == 2:
            raise OSError("injected publish failure")
        return real_replace(source, destination)

    monkeypatch.setattr(input_contract.os, "replace", fail_second_replace)
    with pytest.raises(OSError, match="injected publish failure"):
        input_contract.replace_contract(str(tmp_path))
    assert len(calls) == 3
    assert state_snapshot(tmp_path) == before


def test_missing_or_drifted_contract_fails_validation(tmp_path):
    add_pair(tmp_path, "Alpha")
    with pytest.raises(InputContractError, match="cannot load input contract"):
        input_contract.load_and_validate_contract(str(tmp_path))

    input_contract.create_contract(str(tmp_path))
    write_fastq(
        tmp_path / "Rawdata" / "Alpha" / "Alpha_1.fq.gz",
        read_id="changed",
    )
    with pytest.raises(InputContractError, match="Rawdata differs"):
        input_contract.load_and_validate_contract(str(tmp_path))


def test_header_only_cell_tagging_preserves_fastq_payload(tmp_path):
    input_path = tmp_path / "input.fq.gz"
    output_path = tmp_path / "output.fq.gz"
    code = input_contract.safe_code_for_sample("Alpha_one")
    with gzip.open(input_path, "wb") as handle:
        handle.write(
            b"@instrument_lane_read_ACGTACGT 1:N:0:1\nACGT\n+instrument_lane_read\n____\n"
        )

    assert tag_fastq(str(input_path), str(output_path), code) == 1
    with gzip.open(output_path, "rb") as handle:
        lines = handle.readlines()
    assert lines[0] == (
        "@instrument_lane_read_{}_ACGTACGT 1:N:0:1\n".format(code).encode("ascii")
    )
    assert lines[1:] == [b"ACGT\n", b"+instrument_lane_read\n", b"____\n"]


def test_matrix_round_trip_restores_exact_sample_names_and_order(tmp_path):
    add_pair(tmp_path, "Beta-2")
    add_pair(tmp_path, "Alpha_one")
    contract = input_contract.create_contract(str(tmp_path))
    samples = contract["samples"]
    codes = [sample["safe_code"] for sample in samples]
    matrix = tmp_path / "counts.safe.tsv"
    matrix.write_text(
        "gene\t{}\t{}\nGeneA\t7\t3\n".format(codes[1], codes[0])
    )
    output = tmp_path / "counts.tsv"
    receipt = tmp_path / "counts.cell_contract.tsv"
    rename_count_matrix(
        str(matrix),
        str(tmp_path / "qc" / "input_contract" / "current.json"),
        str(output),
        str(receipt),
    )
    assert output.read_text() == "gene\tAlpha_one\tBeta-2\nGeneA\t3\t7\n"
    assert "cell_column_contract\tPASS" in receipt.read_text()


@pytest.mark.parametrize(
    "header,error",
    [
        ("gene\t{a}\t{a}\n", "duplicated cell codes"),
        ("gene\t{a}\n", "cell-code mismatch"),
        ("gene\t{a}\tUnexpected\n", "cell-code mismatch"),
    ],
)
def test_matrix_rejects_duplicate_missing_and_extra_codes(tmp_path, header, error):
    add_pair(tmp_path, "Alpha")
    add_pair(tmp_path, "Beta")
    contract = input_contract.create_contract(str(tmp_path))
    code = contract["samples"][0]["safe_code"]
    matrix = tmp_path / "counts.safe.tsv"
    matrix.write_text(header.format(a=code) + "GeneA\t1\t2\n")
    with pytest.raises(ValueError, match=error):
        rename_count_matrix(
            str(matrix),
            str(tmp_path / "qc" / "input_contract" / "current.json"),
            str(tmp_path / "counts.tsv"),
            str(tmp_path / "receipt.tsv"),
        )


def test_generated_tsv_is_a_receipt_not_a_user_manifest(tmp_path):
    add_pair(tmp_path, "Alpha_one")
    contract = input_contract.create_contract(str(tmp_path))
    tsv = tmp_path / "qc" / "input_contract" / "discovered_cells.tsv"
    header, row = tsv.read_text().splitlines()
    assert header.startswith("sample_name\tsafe_code\tnaming_scheme")
    assert row.split("\t")[:3] == [
        "Alpha_one",
        contract["samples"][0]["safe_code"],
        "underscore_numeric",
    ]
    current = json.loads(
        (tmp_path / "qc" / "input_contract" / "current.json").read_text()
    )
    assert current == contract
