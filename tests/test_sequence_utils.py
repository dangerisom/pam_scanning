"""Unit tests for the pure sequence-manipulation helpers (no BLAST required)."""

from pam_scanning import library as L


def test_reverse_complement_basic():
    assert L.reverseComplement("ATGC") == "GCAT"
    assert L.reverseComplement("AAAA") == "TTTT"


def test_reverse_complement_preserves_case():
    assert L.reverseComplement("atgc") == "gcat"


def test_reverse_complement_is_involution():
    seq = "ATGCGTACCGGATCAGT"
    assert L.reverseComplement(L.reverseComplement(seq)) == seq


def test_complement_no_reversal():
    assert L.complement("ATGC") == "TACG"


def test_count_mismatches():
    assert L.countMismatches("AAAA", "AAAA") == 0
    assert L.countMismatches("AAAA", "AAGA") == 1
    assert L.countMismatches("ACGT", "TGCA") == 4


def test_fasta_wraps_at_60():
    seq = "A" * 130
    wrapped = L.fasta(seq)
    lines = [ln for ln in wrapped.split("\n") if ln]
    assert lines[0] == "A" * 60
    assert lines[1] == "A" * 60
    assert lines[2] == "A" * 10
    # No wrapped line exceeds 60 characters.
    assert all(len(ln) <= 60 for ln in lines)


def test_mark_silencers_lowercases_unchanged_bases():
    # Identical sequences -> everything lowercased (no mutations highlighted).
    assert L.markSilencers("ATGC", "ATGC") == "atgc"
    # A single change is kept uppercase to flag the silencing mutation.
    assert L.markSilencers("ATTC", "ATGC") == "atTc"
