# How PAM-scanning works: the calculation pipeline

This document explains what the tool computes and why, from the input ORF to the
final primer order. It is the conceptual companion to [`usage.md`](usage.md)
(which lists every flag) and [`blast_setup.md`](blast_setup.md) (BLAST setup).

The whole pipeline runs in one function, `pam_scanning.chimeras.pamscan()`, which
the CLI (`pam-scan`), the GUI (`pam-scan-gui`), and any Python caller share. The
core algorithm lives in `pam_scanning.library`.

## The goal

"PAM scanning" designs, for the codons of an ORF — **every codon by default, or a
chosen subset** — a CRISPR/Cas9 experiment that inserts a sequence (a chimera payload)
at that codon. For each insertion site the tool must produce:

1. a **guide RNA** that directs Cas9 to cut near the insertion codon,
2. a **silenced** version of that guide's target so the edited allele can't be
   re-cut, and
3. the **primers** needed to build the guide plasmid and amplify the payload with
   homology arms for repair.

A guide is only usable if it (a) cuts close enough to the insertion point, (b) can
be silenced without changing the protein, and (c) does not risk cutting elsewhere
in the genome. The pipeline is the sequence of filters that enforces those rules.

> **PAM scanning is always performed in yeast.** The ORF is ported into yeast from
> whatever organism it comes from, so the off-target genome is always a yeast
> genome (default: *S. cerevisiae* BY4741). The gene's source organism never enters
> the off-target step.

## Inputs and the "ORF-plus" sequence

The scan needs the ORF (ATG→stop) plus **100 bp of genomic context on each side**:

- **5′ flank** — 100 bp immediately upstream of the ATG (the `−` side),
- **3′ flank** — 100 bp immediately downstream of the stop (the `+` side).

The driver assembles these into one working sequence:

```
orfPlusSequence = flank5  +  ORF  +  flank3
```

The flanks matter because a guide that edits the **first** or **last** codons must
be able to sit partly outside the ORF. Scanning `orfPlusSequence` (rather than the
bare ORF) is what lets the method reach positions at the very ends of the gene.
The ORF's offset inside this sequence (`orfStartIndex`, always `len(flank5)` = 100)
anchors every codon and cut-site coordinate that follows.

Which codons are actually targeted is then narrowed:

- a **codon selection** restricts insertion to specific residues — as an `.xlsx`,
  a `--codon-positions` list, or the GUI **Pick codons…** picker — or
- a **codon-sampling gap** inserts at every *N*th codon (default 1 = exhaustive).

## Stage 1 — Find PAM sites and candidate guides

`findPamSites()` slides across `orfPlusSequence` and records every SpCas9 PAM on
both strands:

- **Forward strand (`NGG`):** a 23-mer guide = 20 bp protospacer + the `NGG` PAM;
  Cas9 cuts 3 bp 5′ of the PAM.
- **Reverse strand (`CCN`)** (i.e. `NGG` read on the bottom strand): the mirror
  case, with the cut site computed on the reverse strand.

Each guide is keyed by its cut-site coordinate, so downstream stages can ask "how
far is this guide's cut from the codon I want to edit?"

## Stage 2 — Silence the PAM (prevent re-cutting)

After a successful edit, the guide would happily cut the repaired allele again. To
stop that, `tryToPamSilence()` introduces **synonymous (silent) codon
substitutions** that break the PAM (`NGG`/`CCN`) while leaving the encoded protein
unchanged. It first classifies each PAM as **in-frame** or **out-of-frame** with
the ORF reading frame, because that determines which codon(s) overlap the PAM and
which synonymous swaps are available. The bundled yeast codon table (or a
user-supplied one) provides the legal synonymous alternatives.

## Stage 3 — Guide-silence whatever couldn't be PAM-silenced

Some PAMs can't be silenced synonymously (every codon change there would alter the
protein). For those, `guideSilence()` places silent mutations in the **guide body**
instead — as close to the PAM/seed region as possible (where mismatches most
disrupt Cas9), up to a mismatch budget (default 4). Guides that survive **either**
silencing route move on; those that can't be silenced at all are dropped.

## Stage 4 — Screen for off-target cutting with BLAST+

`blastGuides()` runs each candidate guide through **NCBI BLAST+ (`blastn`)** against
the (yeast) host-genome database, then evaluates every hit in Python:

- guides with a credible second cut site in the genome are discarded (**unsafe**);
- guides that are clean but carry marginal matches are flagged as **PAM inclusions**
  — potential, PAM-adjacent off-target sites — at two stringencies (conservative and
  super-conservative), controlled by `--max-pam-inclusions` and
  `--max-pam-inclusion-length`.

Only **safe** guides continue. (The BLAST database is built automatically from the
genome the first time it's used and cached; see [`blast_setup.md`](blast_setup.md).
The `blastn` call itself is fast on the small yeast database — the cost is the
per-hit genome scan in Python, which is parallelized across CPU cores.)

## Stage 5 — Pick the optimal guide per codon

For each target codon, `getOptiGuides()` scores the surviving guides and chooses one:

1. keep only guides whose cut site is within **`maxPamCutGap / 2`** bp of the
   insertion point (default 60/2 = **30 bp**, the empirical proximity rule for
   efficient HDR);
2. among those, rank by **fewest PAM inclusions**, then by **smallest cut-to-codon
   gap**;
3. take the best. A codon with no qualifying guide is recorded as **unscannable**.

The same routine is run across **all** codons (not just the requested ones) to
compute the ORF's overall PAM-scannability and the per-codon cut-gap map used in the
summary plot — a smaller gap means Cas9 cuts closer to the insertion, i.e. higher
expected editing efficiency.

## Stage 6 — Build primers

For each codon that has an optimal guide, the driver emits:

- a **guide primer** — the (silenced-aware) targeting sequence + a suffix that
  amplifies the CRISPR plasmid;
- **insertion primers** (forward and reverse) — left/right **homology arms** taken
  from the *silenced* ORF around the insertion codon, joined to the chimera-payload
  suffixes, trimmed to `--primer-length` (default 100 bp).

Because the homology arms come from the silenced ORF, the repaired allele carries
the silent mutations and is immune to re-cutting.

## Stage 7 — Lay out the order and write QC

`createPrimerOrder()` arranges all guide and insertion primers into a **96- or
384-well plate** (column-major well order) ready to send to a vendor. Alongside it
the driver writes:

- a **guide-solutions workbook** (per codon: cut site, gap, PAM inclusions, original
  and silenced guides, primers),
- **per-codon silenced-ORF FASTAs** (viewable in SnapGene),
- the **scannable-sequence report** and a colour-blind-safe **scannability grid**
  (one cell per codon, coloured by cut gap; when a specific codon set was requested,
  only those codons are reported — unselected codons are greyed out, and a selected
  but inaccessible codon stays white),
- **warnings** for PAM inclusions and unscannable codons,
- copies of the inputs and the assembled `orfPlusContext.fa` for provenance.

See [`usage.md` § Output](usage.md#output) for the exact file layout. A run also
returns a small result dict (`fraction_scannable`, `output_dir`, plot path, and a
text summary) and prints that summary to the console.

## At a glance

| Stage | Function (`pam_scanning.library`) | Removes / produces |
| --- | --- | --- |
| Assemble context | `chimeras.pamscan` | `flank5 + ORF + flank3` |
| 1. Find PAMs | `findPamSites` | all NGG/CCN guides on both strands |
| 2. PAM-silence | `tryToPamSilence` | synonymous PAM knockouts |
| 3. Guide-silence | `guideSilence` | silences the rest; drops the unsilenceable |
| 4. Off-target BLAST | `blastGuides` | drops unsafe guides; flags PAM inclusions |
| 5. Optimal guide | `getOptiGuides` | one guide per codon (or "unscannable") |
| 6. Primers | `chimeras.pamscan` | guide + insertion primers with homology arms |
| 7. Order + QC | `createPrimerOrder`, `calculateScannableSequence`, `plots` | plate order, workbooks, reports, warnings |

Every stage is a filter: a codon reaches the primer order only if a guide passed
proximity, silencing, **and** off-target screening. That is what makes the output
directly orderable.
