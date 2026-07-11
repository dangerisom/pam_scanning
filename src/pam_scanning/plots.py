"""Publication-quality plot of PAM-scannable positions along an ORF.

Renders a compact horizontal track: every codon of the ORF is drawn, coloured by
whether a chimera insertion at that position is reachable by a validated guide
(PAM-scannable) or falls in a gap. Saved as vector PDF and 300-dpi PNG in the QC
folder; the PNG is also shown in the GUI progress console.

matplotlib is imported lazily (with the non-interactive Agg backend) so importing
this module -- or running the CLI without plotting -- stays light and thread-safe.
"""

import os

# Colour-blind-safe pairing on the blue-yellow axis (safe for red-green CVD):
# teal = scannable (good), amber = gap (needs attention).
_SCANNABLE = "#2c7fb8"
_GAP = "#e8a33d"
_TEXT = "#1f2a36"
_MUTED = "#6b7a8d"


def codon_scannability(scannable_sequence):
    """Return a per-codon list of booleans (True = PAM-scannable) from the masked ORF.

    ``scannable_sequence`` is the ORF with unscannable positions blanked to spaces
    (as produced by :func:`pam_scanning.library.calculateScannableSequence`). A
    codon counts as scannable if any of its three nucleotides is unblanked.
    """
    n = len(scannable_sequence) // 3
    return [any(ch != " " for ch in scannable_sequence[3 * c:3 * c + 3]) for c in range(n)]


def _runs(flags, value):
    """Yield (start, length) for each maximal run of *value* in *flags* (0-based)."""
    start = None
    for i, f in enumerate(flags):
        if f == value and start is None:
            start = i
        elif f != value and start is not None:
            yield start, i - start
            start = None
    if start is not None:
        yield start, len(flags) - start


def plot_scannable_positions(qc_path, gene_name, scannable_sequence, fraction,
                             formats=("pdf", "png")):
    """Render the PAM-scannability track for one ORF; return the saved PNG path (or None).

    Writes ``<gene>-scannableMap.<ext>`` into *qc_path* for each requested format.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    flags = codon_scannability(scannable_sequence)
    n = len(flags)
    if n == 0:
        return None
    n_scannable = sum(flags)

    # A tall-enough figure with the track strip in the lower band, leaving a clear
    # header band above for the title, subtitle and legend (fixed positions so they
    # never overlap the bar regardless of ORF length).
    fig = plt.figure(figsize=(10, 2.5))
    ax = fig.add_axes((0.055, 0.34, 0.90, 0.24))   # [left, bottom, width, height]

    # Full ORF as the "gap" baseline, scannable runs drawn over it. Bars are in
    # 1-based residue coordinates: codon c (0-based) spans residues [c+0.5, c+1.5).
    ax.broken_barh([(0.5, n)], (0, 1), facecolors=_GAP, edgecolor="none")
    ax.broken_barh([(s + 0.5, w) for s, w in _runs(flags, True)], (0, 1),
                   facecolors=_SCANNABLE, edgecolor="none")

    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Residue position", fontsize=11)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)

    fig.text(0.055, 0.90, "PAM-scannable positions — %s" % gene_name,
             fontsize=13, fontweight="bold", color=_TEXT, va="top")
    fig.text(0.055, 0.74, "%d of %d codons scannable  (%.1f%% of the ORF)"
             % (n_scannable, n, 100.0 * fraction), fontsize=10.5, color=_MUTED, va="top")
    fig.legend(handles=[Patch(facecolor=_SCANNABLE, label="Scannable"),
                        Patch(facecolor=_GAP, label="Gap (unscannable)")],
               loc="upper right", bbox_to_anchor=(0.955, 0.99), ncol=1, frameon=False,
               fontsize=9.5, handlelength=1.2, labelspacing=0.4)

    png_path = None
    for ext in formats:
        path = os.path.join(qc_path, "%s-scannableMap.%s" % (gene_name, ext))
        fig.savefig(path, dpi=300)
        if ext == "png":
            png_path = path
    plt.close(fig)
    return png_path
