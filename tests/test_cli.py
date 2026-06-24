"""Unit tests for the command-line front-end (no BLAST required).

These cover the new flank-based inputs and the multi-ORF manifest, exercising
the argument/config/manifest plumbing without invoking the BLAST pipeline.
"""

import os

import pytest

from pam_scanning import cli


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(REPO_ROOT, "examples", "fasta")
ORF = os.path.join(FASTA, "S288C_YBL016W_FUS3_coding.fa")
FLANK5 = os.path.join(FASTA, "S288C_YBL016W_FUS3_flank5.fa")
FLANK3 = os.path.join(FASTA, "S288C_YBL016W_FUS3_flank3.fa")


def test_single_orf_flags_build_kwargs():
    base, args = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", "/tmp/genome.fsa", "--gene-name", "Fus3",
    ])
    assert args.manifest is None
    assert base["orf_file_path"] == ORF
    assert base["flank5_file_path"] == FLANK5
    assert base["flank3_file_path"] == FLANK3
    assert base["geneName"] == "Fus3"
    # Untouched defaults survive.
    assert base["primerLength"] == 100
    assert base["localBlastDb"] == "yeast"


def test_validate_requires_both_flanks():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5,
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    with pytest.raises(SystemExit):
        cli._validate(base)  # missing --flank3


def test_validate_passes_for_complete_single_orf():
    base, _ = cli.build_kwargs([
        "--orf", ORF, "--flank5", FLANK5, "--flank3", FLANK3,
        "--genome", ORF, "--gene-name", "Fus3",
    ])
    cli._validate(base)  # should not raise


def test_manifest_parses_rows_and_resolves_paths(tmp_path):
    manifest = tmp_path / "orfs.tsv"
    # Relative paths must resolve against the manifest's own directory.
    (tmp_path / "a.fa").write_text(">a\nATG\n")
    (tmp_path / "u.fa").write_text(">u\nAAA\n")
    (tmp_path / "d.fa").write_text(">d\nTTT\n")
    manifest.write_text(
        "gene\torf\tflank5\tflank3\n"
        "GeneA\ta.fa\tu.fa\td.fa\n"
        "GeneB\ta.fa\tu.fa\td.fa\n"
    )
    orfs = cli._load_manifest(str(manifest))
    assert [o["geneName"] for o in orfs] == ["GeneA", "GeneB"]
    assert orfs[0]["orf_file_path"] == str(tmp_path / "a.fa")
    assert orfs[0]["flank5_file_path"] == str(tmp_path / "u.fa")
    assert orfs[0]["flank3_file_path"] == str(tmp_path / "d.fa")


def test_manifest_accepts_column_aliases(tmp_path):
    manifest = tmp_path / "orfs.tsv"
    manifest.write_text(
        "name\torf\tupstream\tdownstream\tcodon_selection\n"
        "G1\to.fa\tup.fa\tdn.fa\tsel.xlsx\n"
    )
    orfs = cli._load_manifest(str(manifest))
    assert orfs[0]["geneName"] == "G1"
    assert orfs[0]["flank5_file_path"].endswith("up.fa")
    assert orfs[0]["flank3_file_path"].endswith("dn.fa")
    assert orfs[0]["codon_selection_file_path"].endswith("sel.xlsx")


def test_bundled_manifest_example_loads():
    example = os.path.join(REPO_ROOT, "examples", "manifest.tsv")
    orfs = cli._load_manifest(example)
    assert len(orfs) == 1
    assert orfs[0]["geneName"] == "Fus3"
    assert os.path.isfile(orfs[0]["orf_file_path"])
    assert os.path.isfile(orfs[0]["flank5_file_path"])
    assert os.path.isfile(orfs[0]["flank3_file_path"])
