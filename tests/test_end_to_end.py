"""End-to-end regression test for the full pamscan pipeline.

This test is skipped automatically unless NCBI BLAST+ (``blastn`` and
``makeblastdb``) is on PATH. When BLAST+ is available it builds a tiny local
database from the example flanking sequence, runs :func:`pam_scanning.chimeras.pamscan`
on the Fus3 example ORF, and asserts that the expected result artifacts are
produced. It locks the wiring of the pipeline against future refactors.
"""

import os
import shutil
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(REPO_ROOT, "examples", "fasta")

pytestmark = pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="NCBI BLAST+ (blastn/makeblastdb) not available on PATH",
)


def _write_mini_genome(path):
    """Build a single-chromosome genome FASTA in the >chrN/one-line format the
    pipeline's genome parser expects, from the Fus3 flanking example."""
    flanking = os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_flanking.fa")
    seq = "".join(
        line.strip() for line in open(flanking) if not line.startswith(">")
    ).upper()
    with open(path, "w") as fh:
        fh.write(">chr1\n")
        fh.write(seq + "\n")
    return seq


def test_pamscan_generates_outputs(tmp_path):
    from pam_scanning.chimeras import pamscan

    genome = tmp_path / "mini_genome.fsa"
    _write_mini_genome(str(genome))

    db = tmp_path / "minidb"
    subprocess.run(
        ["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)],
        check=True,
        capture_output=True,
    )

    out_dir = tmp_path / "results"
    out_dir.mkdir()

    pamscan(
        orf_file_path=os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_coding.fa"),
        flank5_file_path=os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_flank5.fa"),
        flank3_file_path=os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_flank3.fa"),
        local_genome_file_path=str(genome),
        codon_table_file_path="No file selected",  # use bundled yeast table
        codon_selection_file_path="No file selected",
        geneName="Fus3",
        localBlastDb=str(db),
        guidePrimerForwardSuffix="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG",
        insertPrimerForwardSuffix="GAAGATGTTGTCTGTTGCTCTATGTCATAT",
        insertPrimerReverseSuffix="CTTCTACAACAGACAACGAGATACAGTATA",
        primerLength=100,
        maxPamCutGap=60,
        codonsSamplingGap=1,
        pamInclusionThreshold=5,
        pamInclusionSequenceThreshold=15,
        outputPath=str(out_dir),
    )

    # A single time-stamped run directory should be created.
    runs = [d for d in out_dir.iterdir() if d.is_dir() and "chimera-insertions" in d.name]
    assert len(runs) == 1, "expected exactly one results directory"
    run = runs[0]

    # Core artifacts exist.
    assert (run / "QC" / "Fus3-guideSolutions.xlsx").is_file()
    assert (run / "ORDER" / "Fus3-primerOrder.xlsx").is_file()
    assert (run / "QC" / "Fus3-scannableSequence.txt").is_file()
    assert (run / "BLAST+").is_dir()
    # At least one per-codon solution FASTA was written.
    assert any((run / "QC" / "solutionsFasta").glob("Fus3-*.fa"))
