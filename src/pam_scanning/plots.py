"""Publication-quality grid of PAM-scannable positions along an ORF.

Every codon is a cell in a numbered grid (rows of ``cols`` residues). Inaccessible
positions -- no validated guide reaches them -- are white. Accessible positions are
coloured on a discrete spectrum by the PAM cut gap of their optimal guide (the
distance between the Cas9 cut and the insertion point); a smaller gap means higher
editing efficiency. Saved to the QC folder as vector PDF and 300-dpi PNG, and shown
in the GUI progress console.

matplotlib is imported lazily (non-interactive Agg backend) so importing this module
-- or running the CLI without plotting -- stays light and thread-safe.
"""

import os

_TEXT = "#1f2a36"
_MUTED = "#6b7a8d"
_GRIDLINE = "#cfd6de"


def plot_scannable_positions(qc_path, gene_name, codon_gaps, max_gap, fraction,
                             formats=("pdf", "png"), cols=25, step=5):
    """Render the PAM-scannability grid for one ORF; return the saved PNG path (or None).

    ``codon_gaps`` is a per-codon list (residue = index + 1) holding the optimal
    guide's PAM cut gap, or ``None`` where the codon is inaccessible. ``max_gap`` is
    the largest possible gap (``maxPamCutGap / 2``), used to scale the colour bins.
    Writes ``<gene>-scannableMap.<ext>`` into *qc_path* for each format.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from matplotlib.patches import Patch, Rectangle

    n = len(codon_gaps)
    if n == 0:
        return None
    rows = -(-n // cols)                       # ceil division
    n_accessible = sum(1 for g in codon_gaps if g is not None)

    # Grid of gap values; NaN (masked) = inaccessible or padding -> drawn white.
    grid = np.full((rows, cols), np.nan)
    for i, g in enumerate(codon_gaps):
        if g is not None:
            grid[divmod(i, cols)] = g
    masked = np.ma.masked_invalid(grid)

    # Discrete colour spectrum over the gap, in `step`-bp bins. A perceptually-uniform,
    # colour-blind-safe map (viridis), sub-ranged to avoid the near-white top so a
    # high-gap cell is never confused with a white (inaccessible) cell.
    top = max(step, int(np.ceil(max_gap / step) * step))
    boundaries = np.arange(0, top + step, step)
    n_bins = len(boundaries) - 1
    cmap = mcolors.ListedColormap(cm.viridis(np.linspace(0.08, 0.90, n_bins)))
    cmap.set_bad("white")
    norm = mcolors.BoundaryNorm(boundaries, cmap.N)

    fig_w = min(13.5, cols * 0.34 + 3.4)
    fig_h = min(15.0, rows * 0.34 + 1.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(top=0.84, bottom=0.09, left=0.07, right=0.86)

    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect="equal")

    # Cell borders via minor gridlines (so white cells still read as grid cells).
    ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
    ax.grid(which="minor", color=_GRIDLINE, linewidth=0.6)
    ax.tick_params(which="minor", length=0)

    # Blank the padding cells past the ORF end in the last row (no border) so they
    # read as empty space, not as white "inaccessible" cells.
    filled_last_row = n - (rows - 1) * cols
    if filled_last_row < cols:
        ax.add_patch(Rectangle((filled_last_row - 0.5, rows - 1 - 0.5),
                               cols - filled_last_row, 1.0,
                               facecolor="white", edgecolor="none", zorder=3))

    ax.set_xticks(range(cols))
    ax.set_xticklabels(range(1, cols + 1), fontsize=7)
    ax.set_yticks(range(rows))
    ax.set_yticklabels([r * cols + 1 for r in range(rows)], fontsize=8)
    ax.set_xlabel("Position within row", fontsize=9)
    ax.set_ylabel("Residue at row start", fontsize=9)
    ax.tick_params(axis="both", which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(0.075, 0.955, "PAM-scannable positions — %s" % gene_name,
             fontsize=14, fontweight="bold", color=_TEXT, va="top")
    fig.text(0.075, 0.905, "%d of %d codons accessible  (%.1f%%)"
             % (n_accessible, n, 100.0 * fraction), fontsize=10.5, color=_MUTED, va="top")

    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, ticks=boundaries)
    cbar.set_label("PAM cut gap (bp)\nlower = higher editing efficiency", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.outline.set_visible(False)

    fig.legend(handles=[Patch(facecolor="white", edgecolor=_GRIDLINE, label="Inaccessible")],
               loc="upper right", bbox_to_anchor=(0.995, 0.985), frameon=False, fontsize=9)

    png_path = None
    for ext in formats:
        path = os.path.join(qc_path, "%s-scannableMap.%s" % (gene_name, ext))
        fig.savefig(path, dpi=300)
        if ext == "png":
            png_path = path
    plt.close(fig)
    return png_path
