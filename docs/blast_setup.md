# Setting up NCBI BLAST+ and a local genome database

PAM-scanning uses NCBI BLAST+ to screen candidate guides against the host genome for
off-target Cas9 cutting. You need (1) the BLAST+ executables and (2) a local nucleotide
database built from your host genome.

## 1. Install BLAST+

The simplest route is conda (this is what `environment.yml` does for you):

```bash
conda install -c bioconda blast
```

**Or let PAM-scanning install it for you.** If you run a scan and `blastn` isn't
found, the **GUI** offers to download BLAST+ automatically, streaming progress to
the console. On the command line, add `--install-blast`:

```bash
pam-scan --orf-dir ./orfs --genome genome.fsa --blast-db yeast --install-blast
```

This downloads the official NCBI BLAST+ binaries for your platform into
`~/.pam_scanning/blast` and puts just that folder on `PATH` for the run — **no
conda required**, and nothing is added to any conda environment or your shell
profile. The download is reused on later runs. (It also installs `makeblastdb`,
which you need for step 2 below.)

Alternatively, download the platform installer from the
[NCBI BLAST+ releases](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
and ensure `blastn` and `makeblastdb` are on your `PATH`:

```bash
blastn -version
makeblastdb -version
```

## 2. Obtain a host-genome FASTA

Provide the genome as a FASTA file. The pipeline's genome parser (used for PAM-wobble
checks) expects each chromosome on a single line with a `>chrN` header, for example:

```
>chr1
ACGT...           # entire chromosome 1 on one line
>chr2
ACGT...
```

For *S. cerevisiae* BY4741, `BY4741_Toronto_2012.fsa` is **bundled with the package** and
used by default — you don't need to supply a genome at all for a standard yeast scan. To
scan against a different yeast species/strain/variant, pass your own genome FASTA with
`--genome` (or **Browse Genome sequence** in the GUI).

## 3. The BLAST database is built for you

You normally **do not build a database yourself**. The first time a genome is used, PAM-scanning
runs `makeblastdb` on it automatically, caches the result in `~/.pam_scanning/blastdb` (keyed to
that genome), and reuses it on later runs. This keeps the database in lock-step with the genome
being scanned, so there is nothing to set up: pick a genome (or accept the bundled default) and run.

**Optional override.** If you already have a prebuilt database (e.g. a large shared one), point
PAM-scanning at it and it will skip the auto-build:

```bash
makeblastdb -in /path/to/genome.fsa -dbtype nucl -out /data/yeast   # if building your own
pam-scan ... --blast-db /data/yeast
```

`--blast-db` accepts a path prefix (e.g. `/data/yeast`), a path to any one member file
(e.g. `/data/yeast.nin` — the prefix is taken automatically), or a bare name resolved via
`$BLASTDB`. In the GUI, **Browse database…** selects an existing database and **Use genome
(auto)** returns to the automatic behavior.

## 4. Run

```bash
pam-scan \
    --orf examples/fasta/S288C_YBL016W_FUS3_coding.fa \
    --flank5 examples/fasta/S288C_YBL016W_FUS3_flank5.fa \
    --flank3 examples/fasta/S288C_YBL016W_FUS3_flank3.fa \
    --genome /path/to/BY4741_Toronto_2012.fsa \
    --blast-db yeast \
    --gene-name Fus3 \
    --output ./results
```
