"""Fetch coding DNA sequences (CDS) for UniProt proteins, ready for PAM-scanning.

UniProt distributes *protein* (amino-acid) FASTA, but PAM-scanning operates on the
*nucleotide* ORF (ATG..stop). This module bridges the two deterministically:

    UniProt accession  --(UniProt REST)-->  curated RefSeq mRNA (NM_/XM_)
    RefSeq mRNA         --(NCBI E-utilities)-->  CDS nucleotide sequence

The RefSeq cross-reference is used rather than one of the many EMBL records so the
choice is reproducible: a reviewed UniProt entry names a single curated RefSeq
mRNA, and NCBI's ``fasta_cds_na`` returns that mRNA's CDS directly. Each CDS is
written as ``<gene>_coding.fa`` -- the filename convention
:func:`pam_scanning.cli.discover_orf_folder` recognizes -- so a folder of UniProt
protein FASTAs becomes a folder of PAM-scannable ORFs.

Console entry point: ``pam-scan-fetch-cds``.

The stop codon is retained (PAM-scanning is given the ORF ATG..stop). Only the
network functions touch the internet; parsing and validation are pure and tested
offline.
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# A UniProt FASTA header looks like '>sp|P60709|ACTB_HUMAN ...' or '>tr|A0A...|...';
# capture the middle accession field. Bare accessions are matched directly.
_HEADER_ACCESSION = re.compile(r">?(?:sp|tr)\|([^|]+)\|")
_BARE_ACCESSION = re.compile(r"^[A-Z0-9]+(?:-\d+)?$", re.IGNORECASE)

_FASTA_EXTS = (".fa", ".fasta", ".fna", ".fas")

# Stop codons, to check a fetched CDS really ends at a stop.
_STOP_CODONS = ("TAA", "TAG", "TGA")

# Standard genetic code, used only to pick the RefSeq isoform whose CDS translates
# to the UniProt canonical protein (disambiguation, not sequence design).
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# --- Pure parsing / validation (no network) --------------------------------

def parse_uniprot_accession(header):
    """Extract the UniProt accession from a FASTA header line, or return None.

    Handles the standard '>db|ACCESSION|ENTRY' UniProt header. Any isoform suffix
    (e.g. '-2') is preserved as UniProt returns it.
    """
    match = _HEADER_ACCESSION.search(header or "")
    return match.group(1) if match else None


def accession_from_fasta(path):
    """Return the UniProt accession named in a protein FASTA's first header."""
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                return parse_uniprot_accession(line)
    return None


def select_refseq_mrna(uniprot_json):
    """From a UniProt entry JSON, return (gene_name, [refseq_mrna_ids]).

    RefSeq nucleotide ids are collected from the entry's cross-references, curated
    ``NM_`` records preferred over predicted ``XM_``. The list is sorted so the
    choice downstream is deterministic; the gene name is the entry's primary gene.
    """
    genes = uniprot_json.get("genes") or []
    gene = None
    if genes:
        gene = (genes[0].get("geneName") or {}).get("value")

    curated, predicted = set(), set()
    for xref in uniprot_json.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "RefSeq":
            continue
        for prop in xref.get("properties", []):
            if prop.get("key") != "NucleotideSequenceId":
                continue
            value = (prop.get("value") or "").strip()
            if value.startswith("NM_"):
                curated.add(value)
            elif value.startswith("XM_"):
                predicted.add(value)
    mrnas = sorted(curated) or sorted(predicted)
    return gene, mrnas


def parse_cds_fasta(text):
    """Parse NCBI ``fasta_cds_na`` output into (header, uppercase_sequence).

    Only the first record is returned; a RefSeq mRNA carries exactly one CDS.
    Raises :class:`ValueError` if the response holds no FASTA record (e.g. an
    E-utilities error string).
    """
    header, seq, seen = None, [], False
    for line in text.splitlines():
        if line.startswith(">"):
            if seen:
                break  # start of a second record; stop
            header, seen = line[1:].strip(), True
        elif seen:
            seq.append(line.strip())
    if not seen:
        raise ValueError("No CDS FASTA record in the E-utilities response.")
    return header, "".join(seq).upper()


def validate_cds(sequence):
    """Return a list of human-readable warnings about a CDS (empty if clean)."""
    warnings = []
    if not sequence:
        return ["the CDS is empty"]
    if not sequence.startswith("ATG"):
        warnings.append("does not start with ATG")
    if sequence[-3:] not in _STOP_CODONS:
        warnings.append("does not end in a stop codon")
    if len(sequence) % 3 != 0:
        warnings.append("length %d is not a multiple of 3" % len(sequence))
    non_dna = sorted(set(sequence) - set("ACGTN"))
    if non_dna:
        warnings.append("contains non-DNA character(s): %s" % ", ".join(non_dna))
    return warnings


def translate(dna):
    """Translate a DNA CDS to its protein, stopping at the first stop codon.

    Unknown codons become 'X'. Used only to match a candidate CDS to the UniProt
    canonical protein, so the exact ORF returned is never guessed at.
    """
    residues = []
    for i in range(0, len(dna) - len(dna) % 3, 3):
        residue = _CODON_TABLE.get(dna[i:i + 3], "X")
        if residue == "*":
            break
        residues.append(residue)
    return "".join(residues)


def uniprot_protein_sequence(uniprot_json):
    """Return the canonical protein sequence from a UniProt entry, or None."""
    return (uniprot_json.get("sequence") or {}).get("value")


def sanitize_gene(name):
    """Reduce a gene label to a safe filename stem."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip()).strip("_")
    return cleaned or "gene"


def wrap_fasta(sequence, width=70):
    """Wrap a sequence to fixed-width FASTA body lines."""
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width)) + "\n"


# --- Network -----------------------------------------------------------------

def _http_get(url, params=None, timeout=30):
    """GET a URL and return decoded text, with a descriptive error on failure."""
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "pam_scanning/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        raise RuntimeError("HTTP %s for %s" % (exc.code, url)) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError("Could not reach %s (%s)" % (url, exc.reason)) from exc


def fetch_uniprot_entry(accession, timeout=30):
    """Fetch a UniProt entry (gene name + RefSeq cross-references) as a dict."""
    text = _http_get(
        UNIPROT_URL.format(acc=urllib.parse.quote(accession)),
        params={"fields": "gene_names,sequence,xref_refseq"}, timeout=timeout,
    )
    return json.loads(text)


def fetch_cds_na(refseq_id, timeout=30, email=None):
    """Fetch the CDS nucleotide FASTA for a RefSeq mRNA via NCBI E-utilities."""
    params = {"db": "nuccore", "id": refseq_id,
              "rettype": "fasta_cds_na", "retmode": "text", "tool": "pam_scanning"}
    if email:
        params["email"] = email
    return _http_get(EFETCH_URL, params=params, timeout=timeout)


# --- Orchestration -----------------------------------------------------------

class CdsResult:
    """One resolved CDS: gene, sequence, provenance, and any validation warnings."""

    def __init__(self, accession, gene, refseq_id, sequence, warnings):
        self.accession = accession
        self.gene = gene
        self.refseq_id = refseq_id
        self.sequence = sequence
        self.warnings = warnings


def _cds_of(refseq_id, timeout, email):
    """Fetch and parse the CDS of one RefSeq mRNA into its nucleotide sequence."""
    return parse_cds_fasta(fetch_cds_na(refseq_id, timeout=timeout, email=email))[1]


def _resolve_cds(mrnas, protein, timeout, email):
    """Choose one RefSeq mRNA's CDS, matching the UniProt protein when there's a choice.

    With a single candidate it is used directly. With several, each is translated
    and compared to the UniProt canonical protein, so the isoform is chosen by
    identity rather than by an arbitrary sort. Returns (refseq_id, cds, warnings).
    """
    if len(mrnas) == 1:
        return mrnas[0], _cds_of(mrnas[0], timeout, email), []

    if not protein:
        chosen = mrnas[0]
        return chosen, _cds_of(chosen, timeout, email), [
            "multiple RefSeq mRNAs %s; used %s (no UniProt protein to disambiguate)" % (mrnas, chosen)]

    first_id = first_cds = None
    for refseq_id in mrnas:
        cds = _cds_of(refseq_id, timeout, email)
        if first_cds is None:
            first_id, first_cds = refseq_id, cds
        if translate(cds) == protein:
            note = [] if refseq_id == mrnas[0] else \
                ["selected %s of %s by UniProt protein match" % (refseq_id, mrnas)]
            return refseq_id, cds, note
    return first_id, first_cds, [
        "no RefSeq CDS matched the UniProt canonical protein; used %s of %s" % (first_id, mrnas)]


def fetch_cds_for_accession(accession, timeout=30, email=None):
    """Resolve one UniProt accession to a :class:`CdsResult`.

    One UniProt request plus one NCBI request per RefSeq candidate examined (a
    single request for the common unambiguous case).
    """
    entry = fetch_uniprot_entry(accession, timeout=timeout)
    gene, mrnas = select_refseq_mrna(entry)
    if not mrnas:
        raise LookupError("no RefSeq mRNA cross-reference for %s" % accession)
    refseq_id, sequence, warnings = _resolve_cds(
        mrnas, uniprot_protein_sequence(entry), timeout, email)
    warnings.extend(validate_cds(sequence))
    return CdsResult(accession, gene or accession, refseq_id, sequence, warnings)


def write_coding_fasta(out_dir, result):
    """Write a :class:`CdsResult` as ``<gene>_coding.fa``; return the path."""
    gene = sanitize_gene(result.gene)
    path = os.path.join(out_dir, gene + "_coding.fa")
    header = "> %s | CDS from RefSeq %s (UniProt %s) | %d bp" % (
        gene, result.refseq_id, result.accession, len(result.sequence))
    with open(path, "w") as fh:
        fh.write(header + "\n")
        fh.write(wrap_fasta(result.sequence))
    return path


# --- CLI ---------------------------------------------------------------------

def _gather_accessions(args):
    """Resolve the requested accessions from --protein-dir and/or positional args.

    Returns a list of (accession, source_label) pairs, preserving order and
    dropping duplicates so a folder plus explicit accessions compose cleanly.
    """
    pairs = []
    if args.protein_dir:
        for name in sorted(os.listdir(args.protein_dir)):
            if os.path.splitext(name)[1].lower() not in _FASTA_EXTS:
                continue
            path = os.path.join(args.protein_dir, name)
            accession = accession_from_fasta(path)
            if accession:
                pairs.append((accession, name))
            else:
                print("Warning: no UniProt accession in header of %s; skipped." % name,
                      file=sys.stderr)
    for accession in args.accessions:
        pairs.append((accession, accession))

    seen, unique = set(), []
    for accession, source in pairs:
        if accession not in seen:
            seen.add(accession)
            unique.append((accession, source))
    return unique


def _build_parser():
    p = argparse.ArgumentParser(
        prog="pam-scan-fetch-cds",
        description="Fetch coding DNA sequences (CDS) for UniProt proteins and write "
                    "them as '<gene>_coding.fa' for PAM-scanning folder discovery.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("accessions", nargs="*",
                   help="UniProt accessions (e.g. P60709). Optional if --protein-dir is given.")
    p.add_argument("--protein-dir",
                   help="Folder of UniProt protein FASTAs; the accession is read from "
                        "each file's header.")
    p.add_argument("-o", "--output",
                   help="Directory to write '<gene>_coding.fa' into "
                        "(default: --protein-dir, else the current directory).")
    p.add_argument("--email",
                   help="Contact email passed to NCBI E-utilities (courtesy; recommended "
                        "for large batches).")
    p.add_argument("--timeout", type=int, default=30, help="Per-request timeout (seconds).")
    p.add_argument("--delay", type=float, default=0.34,
                   help="Delay between accessions to respect NCBI rate limits.")
    return p


def main(argv=None):
    """Console-script entry point for ``pam-scan-fetch-cds``."""
    args = _build_parser().parse_args(argv)
    if not args.protein_dir and not args.accessions:
        sys.exit("Error: give at least one accession or --protein-dir.")
    if args.protein_dir and not os.path.isdir(args.protein_dir):
        sys.exit("Error: --protein-dir is not a directory: %s" % args.protein_dir)

    out_dir = args.output or args.protein_dir or "."
    os.makedirs(out_dir, exist_ok=True)

    todo = _gather_accessions(args)
    if not todo:
        sys.exit("Error: no UniProt accessions to fetch.")

    total, written, failed = len(todo), 0, 0
    for i, (accession, source) in enumerate(todo, start=1):
        label = accession if source == accession else "%s (%s)" % (accession, source)
        print("[%d/%d] %s" % (i, total, label))
        try:
            result = fetch_cds_for_accession(accession, timeout=args.timeout, email=args.email)
        except (RuntimeError, LookupError, ValueError) as exc:
            print("  Failed: %s" % exc, file=sys.stderr)
            failed += 1
            continue
        path = write_coding_fasta(out_dir, result)
        written += 1
        note = ("  (warnings: %s)" % "; ".join(result.warnings)) if result.warnings else ""
        print("  %s -> %s [%d bp]%s" % (result.refseq_id, os.path.basename(path),
                                        len(result.sequence), note))
        if i < total and args.delay > 0:
            time.sleep(args.delay)

    print("\nDone: %d written, %d failed, into %s" % (written, failed, out_dir))
    return 1 if failed and not written else 0


if __name__ == "__main__":
    sys.exit(main())
