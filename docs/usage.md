# Usage reference

PAM-scanning can be driven from the command line (`pam-scan`), a Tkinter GUI
(`pam-scan-gui`), or programmatically via `pam_scanning.chimeras.pamscan(**kwargs)`.
All three share the same parameters and produce the same output.

## Parameters

| CLI flag | kwarg | Default | Description |
| --- | --- | --- | --- |
| `--orf` | `orf_file_path` | *(required)* | ORF FASTA, ATG → stop. |
| `--orf-plus` | `orf_plus_buffer_file_path` | *(required)* | ORF flanked by ≥100 bp genomic homology each side. |
| `--genome` | `local_genome_file_path` | *(required)* | Host genome FASTA for off-target checks. |
| `--blast-db` | `localBlastDb` | `yeast` | Name/path of the local BLAST+ database. |
| `--gene-name` | `geneName` | *(required)* | Label used in output filenames. |
| `--codon-table` | `codon_table_file_path` | bundled yeast table | Codon-usage table (`.cusp`-style). |
| `--codon-selection` | `codon_selection_file_path` | *(none)* | `.xlsx` of specific residues to target; overrides sampling. |
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
orf_plus_buffer_file_path = "examples/fasta/S288C_YBL016W_FUS3_flanking.fa"
local_genome_file_path = "/path/to/BY4741_Toronto_2012.fsa"
localBlastDb = "yeast"
geneName = "Fus3"
outputPath = "./results"
codonsSamplingGap = 1
```

```bash
pam-scan --config run.toml
```

## Output

Each run creates a time-stamped directory `‹gene›-chimera-insertions-‹YYYY.MM.DD-HH.MM.SS›/`
under `--output`, containing:

| Path | Contents |
| --- | --- |
| `QC/‹gene›-guideSolutions.xlsx` | Per-codon optimal guide, silenced guide, cut gap, PAM inclusions, and insertion primers. |
| `QC/‹gene›-scannableSequence.txt` | Fraction of the ORF that is PAM-scannable, plus the masked sequence and affected codons. |
| `QC/solutionsFasta/‹gene›-‹codon›.fa` | SnapGene-viewable silenced ORF for each insertion site. |
| `QC/` (copies) | The input ORF and ORF+ FASTA files, for provenance. |
| `ORDER/‹gene›-primerOrder.xlsx` | Plate-laid-out guide + insertion primer order (96- or 384-well). |
| `BLAST+/` | Raw and final `blastn` query/result files for off-target review. |
| `WARNINGS/‹gene›-pamInclusionWarnings*.txt` | Guides carrying potential PAM inclusions (conservative + super-conservative). |
| `WARNINGS/‹gene›-unscannableCodons.txt` | Codons for which no acceptable guide could be found. |
