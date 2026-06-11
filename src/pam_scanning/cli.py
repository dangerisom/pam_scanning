"""Command-line interface for PAM-scanning.

This is a thin, scriptable front-end over :func:`pam_scanning.chimeras.pamscan`.
It exposes exactly the parameters the GUI collects, so GUI and CLI runs are
equivalent and reproducible. Parameters may be given as flags and/or loaded from
a ``--config`` file (JSON or TOML); explicit flags override config values.

Examples
--------
Run an exhaustive scan of an ORF against a local yeast BLAST database::

    pam-scan \\
        --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \\
        --orf-plus examples/fasta/S288C_YBL016W_FUS3_flanking.fa \\
        --genome /path/to/BY4741_Toronto_2012.fsa \\
        --blast-db yeast \\
        --gene-name Fus3 \\
        --output ./results

Or drive everything from a config file::

    pam-scan --config run.toml
"""

import argparse
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


def _build_parser():
    p = argparse.ArgumentParser(
        prog="pam-scan",
        description="Design PAM-scanning chimera-insertion guides and primers for an ORF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs / outputs.
    p.add_argument("--config", help="JSON or TOML file with any of the options below.")
    p.add_argument("--orf", dest="orf_file_path", help="ORF FASTA (ATG..stop).")
    p.add_argument("--orf-plus", dest="orf_plus_buffer_file_path",
                   help="ORF FASTA flanked by genome homology (>=100 bp each side).")
    p.add_argument("--genome", dest="local_genome_file_path",
                   help="Host genome FASTA used for off-target checks.")
    p.add_argument("--codon-table", dest="codon_table_file_path",
                   help="Codon-usage table. Defaults to the bundled yeast table.")
    p.add_argument("--codon-selection", dest="codon_selection_file_path",
                   help="Optional .xlsx of specific insertion sites (overrides sampling).")
    p.add_argument("-o", "--output", dest="outputPath", default=".",
                   help="Directory to write the time-stamped results into.")

    # String parameters.
    p.add_argument("--gene-name", dest="geneName", help="Label used in output filenames.")
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
    """Resolve config file + command-line flags into a full pamscan kwargs dict."""
    args = _build_parser().parse_args(argv)

    kwargs = dict(DEFAULTS)
    if args.config:
        kwargs.update(_load_config(args.config))

    # Flags explicitly provided on the command line override config/defaults.
    for key, value in vars(args).items():
        if key == "config" or value is None:
            continue
        kwargs[key] = value

    # File-path defaults: chimeras.pamscan treats these sentinels as "not provided".
    kwargs.setdefault("codon_table_file_path", "No file selected")
    kwargs.setdefault("codon_selection_file_path", "No file selected")
    kwargs.setdefault("outputPath", ".")
    return kwargs


def _validate(kwargs):
    """Fail early, with clear messages, on missing required inputs."""
    required = {
        "orf_file_path": "--orf",
        "orf_plus_buffer_file_path": "--orf-plus",
        "local_genome_file_path": "--genome",
    }
    missing = [flag for key, flag in required.items()
               if not kwargs.get(key) or kwargs.get(key) == "No file selected"]
    if missing:
        sys.exit("Error: missing required option(s): %s" % ", ".join(sorted(missing)))

    if not kwargs.get("geneName"):
        sys.exit("Error: --gene-name is required.")

    for key, flag in (("orf_file_path", "--orf"),
                      ("orf_plus_buffer_file_path", "--orf-plus"),
                      ("local_genome_file_path", "--genome")):
        path = kwargs[key]
        if not os.path.isfile(path):
            sys.exit("Error: file for %s does not exist: %s" % (flag, path))


def main(argv=None):
    """Console-script entry point for ``pam-scan``."""
    kwargs = build_kwargs(argv)
    _validate(kwargs)
    _check_blast()

    from pam_scanning.chimeras import pamscan

    pamscan(**kwargs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
