# Usage reference

PAM-scanning can be driven from the command line (`pam-scan`), a Tkinter GUI
(`pam-scan-gui`), or programmatically via `pam_scanning.chimeras.pamscan(**kwargs)`.
All three share the same parameters and produce the same output.

## Parameters

| CLI flag | kwarg | Default | Description |
| --- | --- | --- | --- |
| `--orf` | `orf_file_path` | *(required, per ORF)* | ORF FASTA, ATG → stop. |
| `--flank5` | `flank5_file_path` | *(required, per ORF)* | FASTA of the 100 bp immediately **upstream** of the ATG (the `-` side). Lets the scan reach positions at the start of the ORF. |
| `--flank3` | `flank3_file_path` | *(required, per ORF)* | FASTA of the 100 bp immediately **downstream** of the stop (the `+` side). Lets the scan reach positions at the end of the ORF. |
| `--manifest` | *(n/a)* | *(none)* | TSV of ORFs (one row each) for batch runs; see [Multiple ORFs](#multiple-orfs). |
| `--genome` | `local_genome_file_path` | *(required)* | Host genome FASTA for off-target checks. |
| `--blast-db` | `localBlastDb` | `yeast` | Local BLAST+ database. A bare name (e.g. `yeast`) is resolved via `$BLASTDB`; a path prefix or a path to any member file (e.g. `/data/yeast.nin`) is accepted and reduced to the prefix. In the GUI this is a **Browse** button — pick any database file and the prefix path is used. |
| `--gene-name` | `geneName` | *(required, per ORF)* | Label used in output filenames. |
| `--codon-table` | `codon_table_file_path` | bundled yeast table | Codon-usage table (`.cusp`-style). |
| `--codon-selection` | `codon_selection_file_path` | *(none, per ORF)* | `.xlsx` of specific residues to target; overrides sampling. |
| `--output` / `-o` | `outputPath` | `.` | Directory to write the time-stamped run into. |
| `--guide-primer-forward-suffix` | `guidePrimerForwardSuffix` | `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAG` | Suffix that amplifies the CRISPR plasmid after the targeting sequence. |
| `--insert-primer-forward-suffix` | `insertPrimerForwardSuffix` | `GAAGATGTTGTCTGTTGCTCTATGTCATAT` | 5′→3′ insertion-primer suffix (chimera payload, forward). |
| `--insert-primer-reverse-suffix` | `insertPrimerReverseSuffix` | `CTTCTACAACAGACAACGAGATACAGTATA` | 3′→5′ insertion-primer suffix (chimera payload, reverse). |
| `--primer-length` | `primerLength` | `100` | Total length (bp) of the chimera-amplification primers. |
| `--max-pam-cut-gap` | `maxPamCutGap` | `60` | Max bp between sequential PAM cut sites (the empirical 30 bp rule × 2). |
| `--codon-sampling-gap` | `codonsSamplingGap` | `1` | Insert at every Nth codon (1 = exhaustive). Ignored when a codon-selection file is given. |
| `--max-pam-inclusions` | `pamInclusionThreshold` | `5` | Max allowed PAM inclusions per guide solution. |
| `--max-pam-inclusion-length` | `pamInclusionSequenceThreshold` | `15` | Min matched length (bp) counted as a PAM inclusion. |

### Config file

Instead of flags you can pass `--config run.toml` (or `.json`). Keys are the **kwarg**
names from the table above. Flags given on the command line override config values.

```toml
# run.toml
orf_file_path = "examples/fasta/S288C_YBL016W_FUS3_coding.fa"
flank5_file_path = "examples/fasta/S288C_YBL016W_FUS3_flank5.fa"
flank3_file_path = "examples/fasta/S288C_YBL016W_FUS3_flank3.fa"
local_genome_file_path = "/path/to/BY4741_Toronto_2012.fsa"
localBlastDb = "yeast"
geneName = "Fus3"
outputPath = "./results"
codonsSamplingGap = 1
```

```bash
pam-scan --config run.toml
```

## Multiple ORFs

Several ORFs can be scanned in one invocation; each produces its own time-stamped
run directory. There are two ways to supply them.

### Manifest (TSV)

List one ORF per row in a tab-separated **manifest** and supply the shared
parameters with flags or `--config`. Recognized columns: `gene`, `orf`, `flank5`,
`flank3` (required) and `codon_selection` (optional); relative paths resolve
against the manifest's directory. See `examples/manifest.tsv`.

```tsv
gene	orf	flank5	flank3	codon_selection
Fus3	fasta/FUS3_coding.fa	fasta/FUS3_flank5.fa	fasta/FUS3_flank3.fa	codon_selection/Fus3.xlsx
Kss1	fasta/KSS1_coding.fa	fasta/KSS1_flank5.fa	fasta/KSS1_flank3.fa
```

```bash
pam-scan --manifest examples/manifest.tsv \
    --genome /path/to/BY4741_Toronto_2012.fsa --blast-db yeast --output ./results
```

### Folder (`--orf-dir`)

Point at a folder of FASTA files named by convention; the gene name is the part
**before** the role suffix. See `examples/orf_folder/`.

| Role | Suffix (any of) | Type |
| --- | --- | --- |
| ORF | `_coding`, `_orf`, `_cds` | FASTA |
| 5′ flank | `_flank5`, `_5flank`, `_upstream` | FASTA |
| 3′ flank | `_flank3`, `_3flank`, `_downstream` | FASTA |
| codon selection (optional) | `_codonSelection`, `_codons` | `.xlsx` |

```bash
pam-scan --orf-dir examples/orf_folder \
    --genome /path/to/BY4741_Toronto_2012.fsa --blast-db yeast --output ./results
```

Files whose suffix matches no role are ignored with a warning.

### Global vs per-ORF flanks

Flanks can be supplied **per ORF** (each ORF row/file has its own 5′/3′ flank) or
**globally** (one 5′/3′ pair applied to every ORF). For global flanks, omit the
per-ORF flank columns/files and pass `--flank5`/`--flank3` once; they fill in for
any ORF that doesn't specify its own.

```bash
# Global flanks shared by every ORF in the folder:
pam-scan --orf-dir examples/orf_folder \
    --flank5 shared_flank5.fa --flank3 shared_flank3.fa \
    --genome genome.fsa --blast-db yeast --output ./results
```

In the GUI, use **+ Add ORF** to queue ORFs one at a time or **Load folder…** to
discover them, and choose **Per-ORF flanks** or **Global flanks** under *Flank inputs*.

## Output

Each run creates a time-stamped directory `‹gene›-chimera-insertions-‹YYYY.MM.DD-HH.MM.SS›/`
under `--output`, containing:

| Path | Contents |
| --- | --- |
| `QC/‹gene›-guideSolutions.xlsx` | Per-codon optimal guide, silenced guide, cut gap, PAM inclusions, and insertion primers. |
| `QC/‹gene›-scannableSequence.txt` | Fraction of the ORF that is PAM-scannable, plus the masked sequence and affected codons. |
| `QC/solutionsFasta/‹gene›-‹codon›.fa` | SnapGene-viewable silenced ORF for each insertion site. |
| `QC/` (copies) | The input ORF and 5′/3′ flank FASTA files, plus the assembled `‹gene›-orfPlusContext.fa`, for provenance. |
| `ORDER/‹gene›-primerOrder.xlsx` | Plate-laid-out guide + insertion primer order (96- or 384-well). |
| `BLAST+/` | Raw and final `blastn` query/result files for off-target review. |
| `WARNINGS/‹gene›-pamInclusionWarnings*.txt` | Guides carrying potential PAM inclusions (conservative + super-conservative). |
| `WARNINGS/‹gene›-unscannableCodons.txt` | Codons for which no acceptable guide could be found. |
