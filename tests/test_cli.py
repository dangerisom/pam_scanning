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


# --- Folder discovery (--orf-dir) ------------------------------------------

def _touch(path, text=">x\nACGT\n"):
    with open(path, "w") as fh:
        fh.write(text)


def test_discover_folder_groups_by_gene_and_role(tmp_path):
    _touch(tmp_path / "FUS3_coding.fa")
    _touch(tmp_path / "FUS3_flank5.fa")
    _touch(tmp_path / "FUS3_flank3.fa")
    (tmp_path / "FUS3_codonSelection.xlsx").write_bytes(b"PK\x03\x04stub")
    _touch(tmp_path / "KSS1_orf.fasta")
    _touch(tmp_path / "KSS1_upstream.fa")
    _touch(tmp_path / "KSS1_downstream.fa")

    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["FUS3", "KSS1"]
    fus3 = orfs[0]
    assert fus3["orf_file_path"].endswith("FUS3_coding.fa")
    assert fus3["flank5_file_path"].endswith("FUS3_flank5.fa")
    assert fus3["flank3_file_path"].endswith("FUS3_flank3.fa")
    assert fus3["codon_selection_file_path"].endswith("FUS3_codonSelection.xlsx")
    # Alias suffixes resolve to the same roles.
    kss1 = orfs[1]
    assert kss1["orf_file_path"].endswith("KSS1_orf.fasta")
    assert kss1["flank5_file_path"].endswith("KSS1_upstream.fa")
    assert kss1["flank3_file_path"].endswith("KSS1_downstream.fa")


def test_discover_folder_global_flank_layout(tmp_path):
    """Only ORF files present: flanks are expected to come from a global flag."""
    _touch(tmp_path / "GeneA_coding.fa")
    _touch(tmp_path / "GeneB_coding.fa")
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["GeneA", "GeneB"]
    assert all("flank5_file_path" not in o for o in orfs)


def test_discover_folder_reports_unrecognized_and_orphans(tmp_path):
    _touch(tmp_path / "GeneA_coding.fa")
    _touch(tmp_path / "mystery.fa")          # no role suffix
    _touch(tmp_path / "GeneC_flank5.fa")     # flank with no ORF -> orphan
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert [o["geneName"] for o in orfs] == ["GeneA"]
    assert set(skipped) == {"mystery.fa", "GeneC_flank5.fa"}


def test_bundled_orf_folder_example_loads():
    folder = os.path.join(REPO_ROOT, "examples", "orf_folder")
    orfs, skipped = cli.discover_orf_folder(folder)
    assert skipped == []
    assert len(orfs) == 1
    o = orfs[0]
    assert o["geneName"] == "FUS3"
    assert os.path.isfile(o["orf_file_path"])
    assert os.path.isfile(o["flank5_file_path"])
    assert os.path.isfile(o["flank3_file_path"])
    assert os.path.isfile(o["codon_selection_file_path"])


def test_manifest_and_orf_dir_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        cli.main(["--manifest", "m.tsv", "--orf-dir", "d", "--genome", "g.fsa"])


# --- BLAST database prefix -------------------------------------------------

def test_blast_db_prefix_strips_member_extensions():
    assert cli.blast_db_prefix("/data/yeast.nin") == "/data/yeast"
    assert cli.blast_db_prefix("/data/yeast.nsq") == "/data/yeast"
    assert cli.blast_db_prefix("/data/prot.psq") == "/data/prot"
    # Multi-volume member file.
    assert cli.blast_db_prefix("/data/big.00.nhr") == "/data/big"


def test_blast_db_prefix_leaves_names_and_unrelated_paths_alone():
    assert cli.blast_db_prefix("yeast") == "yeast"
    assert cli.blast_db_prefix("/data/yeast") == "/data/yeast"
    assert cli.blast_db_prefix("/genomes/genome.fsa") == "/genomes/genome.fsa"
    assert cli.blast_db_prefix("") == ""


def test_build_kwargs_normalizes_blast_db_path():
    base, _ = cli.build_kwargs(["--blast-db", "/data/yeast.nin", "--orf", ORF])
    assert base["localBlastDb"] == "/data/yeast"
