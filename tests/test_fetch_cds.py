"""Unit tests for the UniProt->CDS fetcher (no network).

Only the pure parsing/validation helpers and the file writer are exercised; the
network functions are not called. UniProt JSON and NCBI FASTA payloads are
represented by small synthetic fixtures mirroring the real response shapes.
"""

import os
import urllib.error

import pytest

from pam_scanning import fetch_cds as fc


class _FakeResponse:
    """Minimal stand-in for a urlopen response usable as a context manager."""

    def __init__(self, body):
        self._body = body

    def read(self):
        return self._body.encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _patch_urlopen(monkeypatch, handler):
    """Route fetch_cds's urlopen through *handler* and make backoff sleeps instant."""
    monkeypatch.setattr(fc.time, "sleep", lambda _s: None)
    monkeypatch.setattr(fc.urllib.request, "urlopen", handler)


# --- HTTP retry / backoff --------------------------------------------------

def test_http_get_retries_transient_5xx_then_succeeds(monkeypatch):
    calls = {"n": 0}

    def handler(request, timeout=30):
        calls["n"] += 1
        if calls["n"] < 3:
            raise urllib.error.HTTPError(request.full_url, 503, "busy", None, None)
        return _FakeResponse("OK")

    _patch_urlopen(monkeypatch, handler)
    assert fc._http_get("http://x", retries=3, backoff=0) == "OK"
    assert calls["n"] == 3   # two failures, then success


def test_http_get_does_not_retry_client_errors(monkeypatch):
    calls = {"n": 0}

    def handler(request, timeout=30):
        calls["n"] += 1
        raise urllib.error.HTTPError(request.full_url, 404, "not found", None, None)

    _patch_urlopen(monkeypatch, handler)
    with pytest.raises(RuntimeError, match="HTTP 404"):
        fc._http_get("http://x", retries=3, backoff=0)
    assert calls["n"] == 1   # 4xx is permanent: no retry


def test_http_get_gives_up_on_persistent_5xx_with_clear_message(monkeypatch):
    calls = {"n": 0}

    def handler(request, timeout=30):
        calls["n"] += 1
        raise urllib.error.HTTPError(request.full_url, 500, "err", None, None)

    _patch_urlopen(monkeypatch, handler)
    with pytest.raises(fc.ServiceUnavailable, match="temporarily unavailable"):
        fc._http_get("http://x", retries=2, backoff=0)
    assert calls["n"] == 3   # retries + 1 attempts


def test_service_unavailable_is_a_runtimeerror():
    # Batch loops catch it specifically to abort, but generic handlers still work.
    assert issubclass(fc.ServiceUnavailable, RuntimeError)


def test_http_get_retries_connection_errors(monkeypatch):
    calls = {"n": 0}

    def handler(request, timeout=30):
        calls["n"] += 1
        if calls["n"] == 1:
            raise urllib.error.URLError("connection reset")
        return _FakeResponse("DATA")

    _patch_urlopen(monkeypatch, handler)
    assert fc._http_get("http://x", retries=3, backoff=0) == "DATA"
    assert calls["n"] == 2


# --- Accession parsing -----------------------------------------------------

def test_parse_accession_from_swissprot_header():
    header = ">sp|P60709|ACTB_HUMAN Actin, cytoplasmic 1 OS=Homo sapiens GN=ACTB"
    assert fc.parse_uniprot_accession(header) == "P60709"


def test_parse_accession_from_trembl_header():
    assert fc.parse_uniprot_accession(">tr|A0A0U1RRL7|A0A0U1RRL7_HUMAN foo") == "A0A0U1RRL7"


def test_parse_accession_returns_none_without_match():
    assert fc.parse_uniprot_accession(">random header no accession") is None


def test_accession_from_fasta_reads_first_header(tmp_path):
    f = tmp_path / "ACTB_P60709.fasta"
    f.write_text(">sp|P60709|ACTB_HUMAN Actin\nMDDDIAALVV\n")
    assert fc.accession_from_fasta(str(f)) == "P60709"


# --- DNA vs protein detection ----------------------------------------------

def test_is_dna_sequence():
    assert fc.is_dna_sequence("ATGACGTNNN")
    assert fc.is_dna_sequence("atgacg")          # case-insensitive
    assert not fc.is_dna_sequence("")            # empty is not DNA
    assert not fc.is_dna_sequence("MDDDIAALVV")  # protein


def test_first_fasta_sequence_reads_one_record(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">a\nATG\nAAA\n>b\nGGG\n")
    assert fc.first_fasta_sequence(str(f)) == "ATGAAA"


def test_fasta_is_protein_distinguishes_dna_and_protein(tmp_path):
    dna = tmp_path / "gene_coding.fa"
    dna.write_text(">g\nATGAAATAG\n")
    protein = tmp_path / "gene.fasta"
    protein.write_text(">sp|P1|G_HUMAN x\nMDDDIAALVV\n")
    assert fc.fasta_is_protein(str(protein)) is True
    assert fc.fasta_is_protein(str(dna)) is False


# --- RefSeq selection from a UniProt entry ---------------------------------

def _uniprot_entry(gene, refseqs):
    return {
        "genes": [{"geneName": {"value": gene}}],
        "uniProtKBCrossReferences": [
            {"database": "RefSeq", "id": "NP_x",
             "properties": [{"key": "NucleotideSequenceId", "value": r}]}
            for r in refseqs
        ],
    }


def test_select_refseq_prefers_curated_nm_and_sorts():
    gene, mrnas = fc.select_refseq_mrna(_uniprot_entry("ACTB", ["NM_001101.5", "XM_999.1"]))
    assert gene == "ACTB"
    assert mrnas == ["NM_001101.5"]          # XM_ dropped when an NM_ exists


def test_select_refseq_falls_back_to_predicted_xm():
    _gene, mrnas = fc.select_refseq_mrna(_uniprot_entry("G", ["XM_2.1", "XM_1.1"]))
    assert mrnas == ["XM_1.1", "XM_2.1"]     # sorted, deterministic


def test_select_refseq_none_when_absent():
    gene, mrnas = fc.select_refseq_mrna({"genes": [{"geneName": {"value": "G"}}]})
    assert gene == "G" and mrnas == []


# --- CDS FASTA parsing -----------------------------------------------------

CDS_FASTA = (
    ">lcl|NM_001101.5_cds_NP_001092.1_1 [gene=ACTB] [location=85..1212]\n"
    "ATGGATGATGAT\n"
    "ACGTAA\n"
)


def test_parse_cds_fasta_first_record_only():
    header, seq = fc.parse_cds_fasta(CDS_FASTA + ">second\nGGGG\n")
    assert "NM_001101.5" in header
    assert seq == "ATGGATGATGATACGTAA"      # joined, upper-cased, second record ignored


def test_parse_cds_fasta_rejects_non_fasta():
    with pytest.raises(ValueError):
        fc.parse_cds_fasta("Error: id not found\n")


# --- CDS validation --------------------------------------------------------

def test_validate_clean_cds_has_no_warnings():
    assert fc.validate_cds("ATG" + "AAA" * 3 + "TAA") == []


def test_validate_flags_start_stop_frame_and_bases():
    warnings = fc.validate_cds("CCGAAT")   # no ATG, no stop, len%3==0 here
    assert any("ATG" in w for w in warnings)
    assert any("stop" in w for w in warnings)
    assert fc.validate_cds("ATGA")         # not a multiple of 3
    assert any("non-DNA" in w for w in fc.validate_cds("ATGXTAA"))


# --- Translation + isoform disambiguation ----------------------------------

def test_translate_stops_at_first_stop_codon():
    assert fc.translate("ATGAAATTTTAAGGG") == "MKF"   # ATG AAA TTT | TAA stop
    assert fc.translate("ATG") == "M"
    assert fc.translate("ATGXX") == "M"               # trailing partial codon ignored


def test_resolve_cds_single_candidate_no_network_for_choice(monkeypatch):
    monkeypatch.setattr(fc, "fetch_cds_na", lambda rid, **k: ">x\nATGAAATAA\n")
    refseq_id, cds, warnings = fc._resolve_cds(["NM_1.1"], "MK", 30, None)
    assert refseq_id == "NM_1.1" and cds == "ATGAAATAA" and warnings == []


def test_resolve_cds_picks_isoform_matching_the_protein(monkeypatch):
    # Two isoforms; only the second translates to the UniProt protein "MF".
    cds_by_id = {"NM_1.1": "ATGAAATAA", "NM_2.2": "ATGTTTTAA"}  # -> "MK" / "MF"
    monkeypatch.setattr(fc, "fetch_cds_na", lambda rid, **k: ">x\n%s\n" % cds_by_id[rid])
    refseq_id, cds, warnings = fc._resolve_cds(["NM_1.1", "NM_2.2"], "MF", 30, None)
    assert refseq_id == "NM_2.2"
    assert cds == "ATGTTTTAA"
    assert any("by UniProt protein match" in w for w in warnings)


def test_resolve_cds_warns_when_no_isoform_matches(monkeypatch):
    cds_by_id = {"NM_1.1": "ATGAAATAA", "NM_2.2": "ATGTTTTAA"}
    monkeypatch.setattr(fc, "fetch_cds_na", lambda rid, **k: ">x\n%s\n" % cds_by_id[rid])
    refseq_id, _cds, warnings = fc._resolve_cds(["NM_1.1", "NM_2.2"], "WWWW", 30, None)
    assert refseq_id == "NM_1.1"  # deterministic fallback to the first
    assert any("no RefSeq CDS matched" in w for w in warnings)


def test_uniprot_protein_sequence_extraction():
    assert fc.uniprot_protein_sequence({"sequence": {"value": "MKF"}}) == "MKF"
    assert fc.uniprot_protein_sequence({}) is None


# --- Filename + writer -----------------------------------------------------

def test_sanitize_gene():
    assert fc.sanitize_gene("ACTB") == "ACTB"
    assert fc.sanitize_gene("NF-κB/p65") == "NF-_B_p65"
    assert fc.sanitize_gene("   ") == "gene"


def test_write_coding_fasta_uses_convention_and_is_discoverable(tmp_path):
    from pam_scanning import cli, chimeras

    result = fc.CdsResult("P60709", "ACTB", "NM_001101.5", "ATGAAATAA", [])
    path = fc.write_coding_fasta(str(tmp_path), result)
    assert os.path.basename(path) == "ACTB_coding.fa"

    # The written file round-trips through the pipeline's own FASTA reader...
    assert chimeras._read_fasta_sequence(path) == "ATGAAATAA"
    # ...and folder discovery recognizes it as ACTB's ORF.
    orfs, skipped = cli.discover_orf_folder(str(tmp_path))
    assert skipped == []
    assert [o["geneName"] for o in orfs] == ["ACTB"]
    assert orfs[0]["orf_file_path"].endswith("ACTB_coding.fa")
