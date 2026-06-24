"""Command-line interface for PAM-scanning.

This is a thin, scriptable front-end over :func:`pam_scanning.chimeras.pamscan`.
It exposes exactly the parameters the GUI collects, so GUI and CLI runs are
equivalent and reproducible. Parameters may be given as flags and/or loaded from
a ``--config`` file (JSON or TOML); explicit flags override config values.

Examples
--------
Scan a single ORF against a local yeast BLAST database::

    pam-scan \\
        --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \\
        --flank5 examples/fasta/S288C_YBL016W_FUS3_flank5.fa \\
        --flank3 examples/fasta/S288C_YBL016W_FUS3_flank3.fa \\
        --genome /path/to/BY4741_Toronto_2012.fsa \\
        --blast-db yeast \\
        --gene-name Fus3 \\
        --output ./results

Scan many ORFs from a tab-separated manifest (one ORF per row), with the shared
parameters supplied by flags or ``--config``::

    pam-scan --manifest orfs.tsv --genome genome.fsa --blast-db yeast --output ./results

The manifest has a header row and one row per ORF. Recognized columns:
``gene``, ``orf``, ``flank5``, ``flank3`` (required) and ``codon_selection``
(optional). Relative paths are resolved against the manifest's directory.
"""

import argparse
import csv
import os
import shutil
import sys


# Parameter defaults shared with the GUI. Keys match the kwargs that
# pam_scanning.chimeras.pamscan expects.
DEFAULTS = {
    "geneName": None,
    "localBlastDb": "yeast",
    "guidePrimerForwardSuffix": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG",
    "insertPrimerForwardSuffix": "GAAGATGTTGTCTGTTGCTCTATGTCATAT",
    "insertPrimerReverseSuffix": "CTTCTACAACAGACAACGAGATACAGTATA",
    "primerLength": 100,
    "maxPamCutGap": 60,
    "codonsSamplingGap": 1,
    "pamInclusionThreshold": 5,
    "pamInclusionSequenceThreshold": 15,
}

# Keys that vary per ORF (supplied by single-ORF flags or by manifest rows).
PER_ORF_KEYS = ("geneName", "orf_file_path", "flank5_file_path", "flank3_file_path",
                "codon_selection_file_path")


def _build_parser():
    p = argparse.ArgumentParser(
        prog="pam-scan",
        description="Design PAM-scanning chimera-insertion guides and primers for one or more ORFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs / outputs.
    p.add_argument("--config", help="JSON or TOML file with any of the options below.")
    p.add_argument("--manifest",
                   help="TSV of ORFs (columns: gene, orf, flank5, flank3, [codon_selection]); "
                        "one ORF per row. Shared params come from flags/--config.")
    p.add_argument("--orf", dest="orf_file_path", help="ORF FASTA (ATG..stop). [single ORF]")
    p.add_argument("--flank5", dest="flank5_file_path",
                   help="5' flank FASTA: 100 bp upstream of the ATG (the '-' side). [single ORF]")
    p.add_argument("--flank3", dest="flank3_file_path",
                   help="3' flank FASTA: 100 bp downstream of the stop (the '+' side). [single ORF]")
    p.add_argument("--genome", dest="local_genome_file_path",
                   help="Host genome FASTA used for off-target checks.")
    p.add_argument("--codon-table", dest="codon_table_file_path",
                   help="Codon-usage table. Defaults to the bundled yeast table.")
    p.add_argument("--codon-selection", dest="codon_selection_file_path",
                   help="Optional .xlsx of specific insertion sites (overrides sampling). [single ORF]")
    p.add_argument("-o", "--output", dest="outputPath", default=".",
                   help="Directory to write the time-stamped results into.")

    # String parameters.
    p.add_argument("--gene-name", dest="geneName", help="Label used in output filenames. [single ORF]")
    p.add_argument("--blast-db", dest="localBlastDb",
                   help="Name of the local BLAST+ database.")
    p.add_argument("--guide-primer-forward-suffix", dest="guidePrimerForwardSuffix")
    p.add_argument("--insert-primer-forward-suffix", dest="insertPrimerForwardSuffix")
    p.add_argument("--insert-primer-reverse-suffix", dest="insertPrimerReverseSuffix")

    # Integer parameters.
    p.add_argument("--primer-length", dest="primerLength", type=int)
    p.add_argument("--max-pam-cut-gap", dest="maxPamCutGap", type=int)
    p.add_argument("--codon-sampling-gap", dest="codonsSamplingGap", type=int)
    p.add_argument("--max-pam-inclusions", dest="pamInclusionThreshold", type=int)
    p.add_argument("--max-pam-inclusion-length", dest="pamInclusionSequenceThreshold", type=int)

    return p


def _load_config(path):
    """Load a JSON or TOML config file into a dict of pamscan kwargs."""
    ext = os.path.splitext(path)[1].lower()
    with open(path, "rb") as fh:
        if ext == ".json":
            import json

            return json.load(fh)
        if ext in (".toml", ".tml"):
            try:
                import tomllib  # Python >= 3.11
            except ModuleNotFoundError:  # pragma: no cover
                try:
                    import tomli as tomllib
                except ModuleNotFoundError:
                    raise SystemExit(
                        "Reading TOML config requires Python 3.11+ or the 'tomli' package."
                    )
            return tomllib.load(fh)
    raise SystemExit("Config file must end in .json or .toml: %s" % path)


# Manifest column name -> pamscan kwarg key. Several spellings are accepted.
_MANIFEST_COLUMNS = {
    "gene": "geneName", "genename": "geneName", "name": "geneName",
    "orf": "orf_file_path",
    "flank5": "flank5_file_path", "5flank": "flank5_file_path", "upstream": "flank5_file_path",
    "flank3": "flank3_file_path", "3flank": "flank3_file_path", "downstream": "flank3_file_path",
    "codon_selection": "codon_selection_file_path",
    "codon-selection": "codon_selection_file_path",
    "codonselection": "codon_selection_file_path",
}


def _load_manifest(path):
    """Parse a TSV manifest into a list of per-ORF kwarg dicts.

    Relative file paths in the manifest are resolved against the manifest's own
    directory so a manifest plus its FASTA files can live together.
    """
    base = os.path.dirname(os.path.abspath(path))

    def resolve(value):
        value = (value or "").strip()
        if not value:
            return None
        return value if os.path.isabs(value) else os.path.normpath(os.path.join(base, value))

    file_keys = {"orf_file_path", "flank5_file_path", "flank3_file_path",
                 "codon_selection_file_path"}
    orfs = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Manifest is empty: %s" % path)
        for raw in reader:
            entry = {}
            for col, value in raw.items():
                key = _MANIFEST_COLUMNS.get((col or "").strip().lower())
                if key is None:
                    continue
                value = (value or "").strip()
                if not value:
                    continue
                entry[key] = resolve(value) if key in file_keys else value
            if not entry:
                continue  # blank line
            orfs.append(entry)
    if not orfs:
        raise SystemExit("Manifest has no ORF rows: %s" % path)
    return orfs


def _check_blast():
    """Ensure the external NCBI BLAST+ 'blastn' executable is available."""
    if shutil.which("blastn") is None:
        sys.exit(
            "Error: 'blastn' (NCBI BLAST+) was not found on your PATH.\n"
            "PAM-scanning requires BLAST+ for off-target evaluation.\n"
            "Install it with conda (recommended):  conda install -c bioconda blast\n"
            "or see https://www.ncbi.nlm.nih.gov/books/NBK279690/"
        )


def build_kwargs(argv=None):
    """Resolve config file + command-line flags into a (base kwargs, args) pair.

    For a single ORF the base is the complete kwargs dict; for a ``--manifest``
    run it is the shared base that each ORF row is merged onto.
    """
    args = _build_parser().parse_args(argv)

    kwargs = dict(DEFAULTS)
    if args.config:
        kwargs.update(_load_config(args.config))

    # Flags explicitly provided on the command line override config/defaults.
    for key, value in vars(args).items():
        if key in ("config", "manifest") or value is None:
            continue
        kwargs[key] = value

    # File-path defaults: chimeras.pamscan treats these sentinels as "not provided".
    kwargs.setdefault("codon_table_file_path", "No file selected")
    kwargs.setdefault("codon_selection_file_path", "No file selected")
    kwargs.setdefault("outputPath", ".")
    return kwargs, args


def _validate(kwargs):
    """Fail early, with clear messages, on missing required inputs for one ORF."""
    required = {
        "orf_file_path": "--orf / manifest 'orf'",
        "flank5_file_path": "--flank5 / manifest 'flank5'",
        "flank3_file_path": "--flank3 / manifest 'flank3'",
        "local_genome_file_path": "--genome",
    }
    missing = [flag for key, flag in required.items()
               if not kwargs.get(key) or kwargs.get(key) == "No file selected"]
    if missing:
        sys.exit("Error: missing required input(s): %s" % ", ".join(sorted(missing)))

    if not kwargs.get("geneName"):
        sys.exit("Error: a gene name is required (--gene-name or manifest 'gene').")

    for key, flag in (("orf_file_path", "--orf"),
                      ("flank5_file_path", "--flank5"),
                      ("flank3_file_path", "--flank3"),
                      ("local_genome_file_path", "--genome")):
        path = kwargs[key]
        if not os.path.isfile(path):
            sys.exit("Error: file for %s does not exist: %s" % (flag, path))


def main(argv=None):
    """Console-script entry point for ``pam-scan`` (single ORF or a --manifest batch)."""
    base, args = build_kwargs(argv)
    _check_blast()

    from pam_scanning.chimeras import pamscan

    if args.manifest:
        orfs = _load_manifest(args.manifest)
        total = len(orfs)
        for i, orf in enumerate(orfs, start=1):
            kwargs = dict(base)
            kwargs.update(orf)                       # per-ORF fields override the shared base
            _validate(kwargs)
            print("=== PAM scan %d/%d: %s ===" % (i, total, kwargs.get("geneName")))
            pamscan(**kwargs)
        return 0

    _validate(base)
    pamscan(**base)
    return 0


if __name__ == "__main__":
    sys.exit(main())
