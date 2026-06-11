"""Tkinter GUI front-end for PAM-scanning.

Collects run parameters through a themed form with hover help and hands them to
:func:`pam_scanning.chimeras.pamscan`. Launch with ``pam-scan-gui``.
"""

import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, font as tkfont
from os import getcwd

# ---------------------------------------------------------------------------
# Palette and field definitions
# ---------------------------------------------------------------------------

BG = "#eef2f5"          # window background
CARD = "#ffffff"        # section card background
HEADER_BG = "#1f3a5f"   # header banner
HEADER_FG = "#ffffff"
SUBTITLE_FG = "#c3d0de"
TEXT = "#1f2a36"
MUTED = "#6b7a8d"
ACCENT = "#2e7d32"      # primary action (Run)
ACCENT_ACTIVE = "#256628"
TIP_BG = "#fffbe6"
TIP_BORDER = "#c9b458"
PLACEHOLDER = "No file selected"

# Each tuple: (kwarg key, browse-button label, tooltip)
FILE_FIELDS = [
    ("orf_file_path", "ORF",
     "Open the open reading frame (ORF) FASTA file for the gene of interest. "
     "The ORF sequence should begin with the ATG start codon and end with a stop codon."),
    ("orf_plus_buffer_file_path", "ORF+",
     "Open the ORF FASTA file flanked by at least ~100 bp (ideally up to 1000 bp) of "
     "genomic homology on each side."),
    ("local_genome_file_path", "Genome sequence",
     "The host genome sequence in FASTA format used for off-target evaluation. For a yeast "
     "PAM scan, this is the full yeast genome."),
    ("codon_table_file_path", "Codon table (optional)",
     "Codon-usage table for the host genome. Leave unset to use the bundled yeast table "
     "(yeast_64_1_1_all_nuclear.cusp.txt)."),
    ("codon_selection_file_path", "Codon selection (optional)",
     "An .xlsx file listing specific chimera insertion points. Providing this overrides the "
     "codon sampling frequency parameter below."),
]

# Each tuple: (kwarg key, label, default, tooltip)
STRING_FIELDS = [
    ("geneName", "Gene name", "MFG",
     "A short gene-name label used in the names of the generated output files."),
    ("localBlastDb", "Local BLAST database", "yeast",
     "The name (or path) of your local BLAST+ database. To create one, install BLAST+ "
     "(https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and run makeblastdb."),
    ("guidePrimerForwardSuffix", "Guide primer: forward suffix",
     "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG",
     "Forward primer sequence that amplifies the CRISPR plasmid starting after the unique "
     "targeting sequence (UnTS). It is prepended with each 20 bp guide to build the guide "
     "plasmid via the around-the-world protocol."),
    ("insertPrimerForwardSuffix", "Insert primer: forward suffix",
     "GAAGATGTTGTCTGTTGCTCTATGTCATAT",
     "The 5' to 3' top-strand insertion-primer suffix concatenated downstream of the PAM-scan "
     "chimera-site genome homology. With the reverse suffix below it amplifies the chimeric "
     "insert for each PAM-scan site."),
    ("insertPrimerReverseSuffix", "Insert primer: reverse suffix",
     "CTTCTACAACAGACAACGAGATACAGTATA",
     "The 3' to 5' bottom-strand insertion-primer suffix concatenated upstream of the PAM-scan "
     "chimera-site genome homology. With the forward suffix above it amplifies the chimeric "
     "insert for each PAM-scan site."),
]

# Each tuple: (kwarg key, label, default, tooltip)
INT_FIELDS = [
    ("primerLength", "Primer length (bp)", 100,
     "Total length of the DNA primers used to amplify the chimera insert."),
    ("maxPamCutGap", "Max gap between two PAMs (bp)", 60,
     "Sequential PAM sites separated by more than this distance create PAM-scanning gaps with "
     "lower editing efficiency. The default of 60 reflects the empirical 30 bp rule (chimeric "
     "insertion >30 bp from the Cas9 cut site is less efficient)."),
    ("codonsSamplingGap", "Codon sampling frequency", 1,
     "Insert the chimera at every Nth codon. 1 = exhaustively after every codon, 2 = every "
     "other codon, and so on. Ignored when a codon selection file is provided."),
    ("pamInclusionThreshold", "Max PAM inclusions", 5,
     "Slightly shorter guides ('PAM inclusions') may still enable cutting; these are tracked "
     "across solutions. 5 generates ample guide designs while keeping PAM inclusions minimal."),
    ("pamInclusionSequenceThreshold", "Max PAM inclusion length (bp)", 15,
     "The minimum matched length counted as a PAM inclusion. With the default of 15, candidate "
     "inclusions longer than 15 bp are counted."),
]


# ---------------------------------------------------------------------------
# Hover tooltip
# ---------------------------------------------------------------------------

class Tooltip:
    """A lightweight hover tooltip that pops up next to a widget after a delay."""

    def __init__(self, widget, text, *, delay=350, wraplength=380):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wraplength = wraplength
        self._after_id = None
        self._tip = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after_id is not None:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip or not self.text:
            return
        x = self.widget.winfo_rootx()
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry("+%d+%d" % (x, y))
        self._tip.attributes("-topmost", True)
        frame = tk.Frame(self._tip, background=TIP_BORDER, bd=0)
        frame.pack()
        label = tk.Label(
            frame, text=self.text, justify="left", background=TIP_BG, foreground=TEXT,
            wraplength=self.wraplength, padx=10, pady=7,
        )
        label.pack(padx=1, pady=1)

    def _hide(self, _event=None):
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

def main():
    root = tk.Tk()
    root.title("PAM Scanning")
    root.configure(bg=BG)
    root.minsize(900, 820)

    # Fonts.
    base = tkfont.nametofont("TkDefaultFont")
    family = base.actual("family")
    # Pick a monospace family that exists on this platform; fall back to Tk's
    # built-in fixed font (which Tk maps to a sensible default everywhere).
    available = set(tkfont.families(root))
    mono_family = next(
        (f for f in ("Menlo", "Consolas", "DejaVu Sans Mono", "Courier New")
         if f in available),
        tkfont.nametofont("TkFixedFont").actual("family"),
    )
    f_title = tkfont.Font(family=family, size=20, weight="bold")
    f_subtitle = tkfont.Font(family=family, size=11)
    f_section = tkfont.Font(family=family, size=12, weight="bold")
    f_label = tkfont.Font(family=family, size=11)
    f_mono = tkfont.Font(family=mono_family, size=10)
    f_button = tkfont.Font(family=family, size=11)
    f_run = tkfont.Font(family=family, size=13, weight="bold")

    # Theme + styles.
    style = ttk.Style()
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    style.configure("TFrame", background=BG)
    style.configure("Card.TFrame", background=CARD)
    style.configure("Header.TFrame", background=HEADER_BG)
    style.configure("Title.TLabel", background=HEADER_BG, foreground=HEADER_FG, font=f_title)
    style.configure("Subtitle.TLabel", background=HEADER_BG, foreground=SUBTITLE_FG, font=f_subtitle)
    style.configure("Section.TLabel", background=BG, foreground=HEADER_BG, font=f_section)
    style.configure("Field.TLabel", background=CARD, foreground=TEXT, font=f_label)
    style.configure("Path.TLabel", background=CARD, foreground=MUTED, font=f_label)
    style.configure("Status.TLabel", background=BG, foreground=MUTED, font=f_subtitle)
    style.configure("TEntry", fieldbackground="#ffffff", padding=4)
    style.configure("Browse.TButton", font=f_button, padding=(10, 4))
    style.configure(
        "Accent.TButton", font=f_run, foreground="#ffffff", background=ACCENT,
        padding=(26, 12), borderwidth=0,
    )
    style.map(
        "Accent.TButton",
        background=[("active", ACCENT_ACTIVE), ("disabled", "#9bbf9d")],
        foreground=[("disabled", "#eef2f5")],
    )

    # Integer-entry validation.
    def _validate_int(proposed):
        return proposed == "" or proposed.isdigit()

    vint = (root.register(_validate_int), "%P")

    # State variables.
    file_vars = {key: tk.StringVar(value=PLACEHOLDER) for key, _, _ in FILE_FIELDS}
    string_vars = {key: tk.StringVar(value=default) for key, _, default, _ in STRING_FIELDS}
    int_vars = {key: tk.StringVar(value=str(default)) for key, _, default, _ in INT_FIELDS}
    output_var = tk.StringVar(value=getcwd())
    status_var = tk.StringVar(value="Ready.")

    # --- Scrollable body so the form fits any screen ---------------------
    outer = ttk.Frame(root, style="TFrame")
    outer.pack(fill="both", expand=True)

    canvas = tk.Canvas(outer, background=BG, highlightthickness=0)
    scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
    body = ttk.Frame(canvas, style="TFrame")
    body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    body_window = canvas.create_window((0, 0), window=body, anchor="nw")
    canvas.bind("<Configure>", lambda e: canvas.itemconfigure(body_window, width=e.width))
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta or (120 if event.num == 4 else -120)) / 120), "units")

    canvas.bind_all("<MouseWheel>", _on_mousewheel)
    canvas.bind_all("<Button-4>", _on_mousewheel)
    canvas.bind_all("<Button-5>", _on_mousewheel)

    # --- Header banner ---------------------------------------------------
    header = ttk.Frame(body, style="Header.TFrame")
    header.pack(fill="x")
    inner_h = ttk.Frame(header, style="Header.TFrame")
    inner_h.pack(fill="x", padx=28, pady=(20, 18))
    ttk.Label(inner_h, text="PAM Scanning", style="Title.TLabel").pack(anchor="w")
    ttk.Label(
        inner_h,
        text="Design CRISPR guides and primers for ORF-wide chimera insertion.",
        style="Subtitle.TLabel",
    ).pack(anchor="w", pady=(4, 0))

    content = ttk.Frame(body, style="TFrame")
    content.pack(fill="both", expand=True, padx=24, pady=18)

    def attach_tip(anchor, text, *widgets):
        """Reveal `text` as a hover tooltip over the anchor and any extra row widgets."""
        tip = Tooltip(anchor, text)
        for w in widgets:
            w.bind("<Enter>", tip._schedule, add="+")
            w.bind("<Leave>", tip._hide, add="+")
            w.bind("<ButtonPress>", tip._hide, add="+")
        return tip

    def section(title, weighted_col):
        ttk.Label(content, text=title, style="Section.TLabel").pack(anchor="w", pady=(14, 6))
        card = ttk.Frame(content, style="Card.TFrame")
        card.pack(fill="x")
        card.columnconfigure(weighted_col, weight=1)
        return card

    def add_file_row(card, row, key, button_label, tip, var):
        path_lbl = ttk.Label(card, textvariable=var, style="Path.TLabel",
                             wraplength=560, anchor="w", justify="left")
        path_lbl.grid(row=row, column=0, sticky="we", padx=(14, 6), pady=10)

        def browse(v=var):
            chosen = filedialog.askopenfilename(initialdir=getcwd())
            if chosen:
                v.set(chosen)

        btn = ttk.Button(card, text="Browse  " + button_label, style="Browse.TButton",
                         command=browse, width=26)
        btn.grid(row=row, column=1, sticky="e", padx=(6, 14), pady=10)
        attach_tip(path_lbl, tip, btn)

    def add_entry_row(card, row, key, label, tip, var, *, mono=False, validate=False):
        lbl = ttk.Label(card, text=label, style="Field.TLabel")
        lbl.grid(row=row, column=0, sticky="w", padx=(14, 6), pady=8)
        kwargs = {"textvariable": var, "font": f_mono if mono else f_label}
        if validate:
            kwargs["validate"] = "key"
            kwargs["validatecommand"] = vint
            kwargs["justify"] = "right"
            kwargs["width"] = 12
        entry = ttk.Entry(card, **kwargs)
        if validate:
            entry.grid(row=row, column=1, sticky="e", padx=(6, 14), pady=8)
        else:
            entry.grid(row=row, column=1, sticky="we", padx=(6, 14), pady=8)
        attach_tip(lbl, tip, entry)
        return entry

    # --- Input files -----------------------------------------------------
    card = section("Input files", 0)
    for r, (key, blabel, tip) in enumerate(FILE_FIELDS):
        add_file_row(card, r, key, blabel, tip, file_vars[key])

    # --- Sequence & primer settings -------------------------------------
    card = section("Sequence & primer settings", 1)
    for r, (key, label, _default, tip) in enumerate(STRING_FIELDS):
        add_entry_row(card, r, key, label, tip, string_vars[key],
                      mono=(key != "geneName" and key != "localBlastDb"))

    # --- Scan parameters -------------------------------------------------
    card = section("Scan parameters", 1)
    for r, (key, label, _default, tip) in enumerate(INT_FIELDS):
        add_entry_row(card, r, key, label, tip, int_vars[key], validate=True)

    # --- Output ----------------------------------------------------------
    card = section("Output", 0)
    out_lbl = ttk.Label(card, textvariable=output_var, style="Path.TLabel",
                        wraplength=560, anchor="w", justify="left")
    out_lbl.grid(row=0, column=0, sticky="we", padx=(14, 6), pady=10)

    def browse_dir():
        chosen = filedialog.askdirectory(initialdir=getcwd())
        if chosen:
            output_var.set(chosen)

    out_btn = ttk.Button(card, text="Choose directory", style="Browse.TButton",
                         command=browse_dir, width=26)
    out_btn.grid(row=0, column=1, sticky="e", padx=(6, 14), pady=10)
    attach_tip(out_lbl, "Directory where the time-stamped PAM-scan output folder is written.", out_btn)

    # --- Action bar ------------------------------------------------------
    action = ttk.Frame(content, style="TFrame")
    action.pack(fill="x", pady=(22, 6))
    run_button = ttk.Button(action, text="▶  Run PAM Scan", style="Accent.TButton")
    run_button.pack()
    ttk.Label(content, textvariable=status_var, style="Status.TLabel").pack(anchor="center", pady=(8, 4))

    # --- Run logic -------------------------------------------------------
    def collect_kwargs():
        kwargs = {key: var.get() for key, var in file_vars.items()}
        kwargs.update({key: var.get() for key, var in string_vars.items()})
        for key, var in int_vars.items():
            text = var.get().strip()
            if not text:
                messagebox.showerror("Invalid input", "Please enter a value for '%s'." % key)
                return None
            kwargs[key] = int(text)
        out = output_var.get().strip()
        kwargs["outputPath"] = out if out and out != PLACEHOLDER else "."
        return kwargs

    def validate_required(kwargs):
        required = {
            "orf_file_path": "ORF",
            "orf_plus_buffer_file_path": "ORF+",
            "local_genome_file_path": "Genome sequence",
        }
        missing = [name for key, name in required.items()
                   if not kwargs[key] or kwargs[key] == PLACEHOLDER]
        if missing:
            messagebox.showerror("Missing input", "Please select: " + ", ".join(missing) + ".")
            return False
        if not kwargs["geneName"].strip():
            messagebox.showerror("Missing input", "Please enter a gene name.")
            return False
        return True

    def finish(error, out_path):
        run_button.config(state="normal")
        if error is not None:
            status_var.set("Error — see message.")
            messagebox.showerror("PAM scan failed", str(error))
        else:
            status_var.set("Done. Output written under: " + out_path)
            messagebox.showinfo("PAM scan complete", "Output written under:\n" + out_path)

    def run_scan():
        kwargs = collect_kwargs()
        if kwargs is None or not validate_required(kwargs):
            return
        run_button.config(state="disabled")
        status_var.set("Running PAM scan… progress is printed to the console.")
        out_path = kwargs["outputPath"]

        def worker():
            try:
                from pam_scanning.chimeras import pamscan
                pamscan(**kwargs)
                root.after(0, lambda: finish(None, out_path))
            except Exception as exc:  # surface any failure back on the UI thread
                root.after(0, lambda e=exc: finish(e, out_path))

        threading.Thread(target=worker, daemon=True).start()

    run_button.config(command=run_scan)

    root.mainloop()


if __name__ == "__main__":
    main()
