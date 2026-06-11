# PAM-scanning

Design CRISPR/Cas9 guide RNAs and chimera-insertion primers for **every codon of an
ORF** ("PAM scanning"). Given an open reading frame and its genomic context, the tool:

1. finds all NGG/CCN PAM sites across the ORF and the guides that target them,
2. **silences** each PAM site (or, failing that, the guide body) with synonymous
   codon substitutions so re-cutting is prevented after editing,
3. uses **NCBI BLAST+** to screen guides for off-target potential against the host
   genome (including conservative "PAM-inclusion" tracking),
4. selects the optimal guide for each insertion codon, and
5. emits ready-to-order guide and insertion primers (96-/384-well plate layout),
   QC FASTA files, a scannable-sequence report, and off-target warnings.

This is the reference implementation accompanying the PAM-scanning manuscript.

## Installation

PAM-scanning depends on the external **NCBI BLAST+** toolkit (`blastn`), which is most
easily obtained through conda. The recommended path is therefore conda.

### conda (recommended)

```bash
git clone https://github.com/USERNAME/pam-scanning.git
cd pam-scanning
conda env create -f environment.yml
conda activate pam-scanning
```

This installs Python, `openpyxl`, **BLAST+**, and the `pam-scanning` package itself
(editable). Once a [Bioconda](https://bioconda.github.io/) release is published you will
also be able to run:

```bash
conda install -c bioconda -c conda-forge pam-scanning
```

### pip

```bash
pip install pam-scanning
```

`pip` does **not** install BLAST+ — install it separately (e.g.
`conda install -c bioconda blast`, or from the
[NCBI BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)).

## Prerequisites: a local BLAST database

You need a local BLAST+ nucleotide database for your host genome. See
[`docs/blast_setup.md`](docs/blast_setup.md) for step-by-step instructions
(`makeblastdb`, where the genome FASTA goes, and how the database name maps to the
`--blast-db` option).

## Usage

### Command line

```bash
pam-scan \
    --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \
    --orf-plus examples/fasta/S288C_YBL016W_FUS3_flanking.fa \
    --genome /path/to/BY4741_Toronto_2012.fsa \
    --blast-db yeast \
    --gene-name Fus3 \
    --output ./results
```

All parameters can also be supplied via `--config run.toml` (or `run.json`); explicit
flags override config values. Run `pam-scan --help` for the full list. The yeast codon
table is bundled and used by default; supply `--codon-table` to override it, and
`--codon-selection` to restrict insertion to specific residues. See
[`docs/usage.md`](docs/usage.md) for the full parameter reference and a description of
the output directory.

### Graphical interface

```bash
pam-scan-gui
```

A Tkinter form that collects the same parameters and runs the identical pipeline.

## Inputs

| Input | Description |
| --- | --- |
| ORF FASTA (`--orf`) | The coding sequence, ATG → stop. |
| ORF+ FASTA (`--orf-plus`) | The same ORF flanked by ≥100 bp of genomic sequence on each side. |
| Genome FASTA (`--genome`) | Host genome used for off-target evaluation (also the source for your BLAST DB). |
| Codon table | Codon-usage table; the bundled yeast table is used if omitted. |
| Codon selection (optional) | `.xlsx` listing specific residues to target. |

## Development

```bash
pip install -e ".[test]"
pytest
```

The unit tests run without BLAST+. The end-to-end regression test
(`tests/test_end_to_end.py`) is skipped automatically unless `blastn`/`makeblastdb`
are on your `PATH`.

## Citation

If you use this software, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE) © 2026 Daniel Isom (Isom Lab)
