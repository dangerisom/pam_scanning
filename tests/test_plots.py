"""Tests for the scannability summary and grid plot (matplotlib runs headless via Agg)."""

from pam_scanning import chimeras, plots


def test_plot_grid_writes_pdf_and_png(tmp_path):
    # 120 codons; a few inaccessible (None), the rest with a PAM cut gap in 0..29.
    gaps = [None if c in (5, 6, 7, 40, 41) else (c % 30) for c in range(120)]
    png = plots.plot_scannable_positions(str(tmp_path) + "/", "GeneX", gaps, 30, 115 / 120)
    assert png == str(tmp_path / "GeneX-scannableMap.png")
    assert (tmp_path / "GeneX-scannableMap.png").is_file()
    assert (tmp_path / "GeneX-scannableMap.pdf").is_file()


def test_plot_empty_orf_returns_none(tmp_path):
    assert plots.plot_scannable_positions(str(tmp_path) + "/", "GeneX", [], 30, 0.0) is None


def test_plot_all_inaccessible_still_renders(tmp_path):
    # Every codon white; must not crash (no colour data).
    png = plots.plot_scannable_positions(str(tmp_path) + "/", "GeneX", [None] * 40, 30, 0.0)
    assert (tmp_path / "GeneX-scannableMap.png").is_file()
    assert png is not None


def test_build_summary_contains_key_fields():
    s = chimeras._build_summary("Fus3", 354, 1062, 1.0, 354, 354, 0, 95, 708, 0, "/out/dir")
    assert "PAM-scan summary: Fus3" in s
    assert "354 codons (1062 bp)" in s
    assert "100.0% of the ORF" in s
    assert "354 of 354 requested" in s
    assert "/out/dir" in s
