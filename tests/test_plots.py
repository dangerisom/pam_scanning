"""Tests for the scannability summary and plot (matplotlib runs headless via Agg)."""

from pam_scanning import chimeras, plots


def test_codon_scannability_from_masked_sequence():
    # codon 0 = "ACG" (scannable), codon 1 = "   " (gap), codon 2 = "ACG" (scannable)
    assert plots.codon_scannability("ACG" + "   " + "ACG") == [True, False, True]


def test_runs_finds_maximal_runs():
    assert list(plots._runs([True, True, False, True], True)) == [(0, 2), (3, 1)]
    assert list(plots._runs([False, False], True)) == []
    assert list(plots._runs([True, True, True], True)) == [(0, 3)]


def test_plot_scannable_positions_writes_pdf_and_png(tmp_path):
    seq = list("ACG" * 20)                    # 20 codons
    for c in (5, 6, 7):                        # a gap
        for k in range(3):
            seq[3 * c + k] = " "
    png = plots.plot_scannable_positions(str(tmp_path) + "/", "GeneX", "".join(seq), 0.85)
    assert png == str(tmp_path / "GeneX-scannableMap.png")
    assert (tmp_path / "GeneX-scannableMap.png").is_file()
    assert (tmp_path / "GeneX-scannableMap.pdf").is_file()


def test_plot_empty_orf_returns_none(tmp_path):
    assert plots.plot_scannable_positions(str(tmp_path) + "/", "GeneX", "", 0.0) is None


def test_build_summary_contains_key_fields():
    s = chimeras._build_summary("Fus3", 354, 1062, 1.0, 354, 354, 0, 95, 708, 0, "/out/dir")
    assert "PAM-scan summary: Fus3" in s
    assert "354 codons (1062 bp)" in s
    assert "100.0% of the ORF" in s
    assert "354 of 354 requested" in s
    assert "/out/dir" in s
