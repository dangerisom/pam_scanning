# Setting up NCBI BLAST+ and a local genome database

PAM-scanning uses NCBI BLAST+ to screen candidate guides against the host genome for
off-target Cas9 cutting. You need (1) the BLAST+ executables and (2) a local nucleotide
database built from your host genome.

## 1. Install BLAST+

The simplest route is conda (this is what `environment.yml` does for you):

```bash
conda install -c bioconda blast
```

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

For *S. cerevisiae* BY4741 we use `BY4741_Toronto_2012.fsa`. The genome file is **not**
shipped with this repository — download it from your preferred source (e.g.
[SGD](https://www.yeastgenome.org/)) and keep it somewhere stable; you pass its path with
`--genome`.

## 3. Build the local BLAST database

```bash
makeblastdb -in /path/to/BY4741_Toronto_2012.fsa -dbtype nucl -out yeast
```

This produces a set of `yeast.*` index files. The name you give `-out` (here `yeast`) is
exactly what you pass to PAM-scanning as `--blast-db` (CLI) or "Local BLAST database"
(GUI). Keep the database files together; point `--blast-db` at the full path or run from
the directory that contains them.

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
