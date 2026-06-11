"""Smoke/structure tests for PAM-site discovery and silencing.

These exercise the silencing machinery on the real, properly-buffered Fus3 example
ORF without requiring BLAST. They assert structural invariants (guides are the
right length and orientation, silenced guides carry silent mutations, the
partition is clean) rather than hand-derived biology.
"""

import os

import pytest

from pam_scanning import library as L
from pam_scanning.chimeras import default_codon_table_path

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(REPO_ROOT, "examples", "fasta")


def _read_fasta(path):
    return "".join(
        line.strip() for line in open(path) if not line.startswith(">")
    ).upper()


@pytest.fixture(scope="module", autouse=True)
def _codon_table():
    # tryToPamSilence / guideSilence rely on the module-level codon tables.
    L.set_codon_table(default_codon_table_path())


@pytest.fixture(scope="module")
def fus3():
    orf = _read_fasta(os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_coding.fa"))
    orf_plus = _read_fasta(os.path.join(EXAMPLES, "S288C_YBL016W_FUS3_flanking.fa"))
    assert orf in orf_plus, "example ORF must be a substring of its flanking sequence"
    start = orf_plus.index(orf)
    end = start + len(orf)
    return orf, orf_plus, start, end


def test_find_pam_sites_orientation(fus3):
    _orf, orf_plus, start, end = fus3
    forward, reverse = L.findPamSites(orf_plus, start, end)
    assert forward and reverse, "Fus3 should contain both forward and reverse PAM sites"
    # Forward guides are 23 bp ending in GG; reverse guides are 23 bp starting with CC.
    assert all(g.endswith("GG") and len(g) == 23 for g in forward.values())
    assert all(g.startswith("CC") and len(g) == 23 for g in reverse.values())


def test_silencing_keeps_guide_length_and_flags_mutations(fus3):
    _orf, orf_plus, start, end = fus3
    forward, reverse = L.findPamSites(orf_plus, start, end)
    guides = dict(forward)
    guides.update(reverse)

    silenced, unsilenced = L.tryToPamSilence(orf_plus, start, end, guides)
    assert silenced, "expected at least some PAM-silenceable guides in Fus3"

    for key, (silenced_guide, _silenced_orf) in silenced.items():
        original = guides[key]
        assert len(silenced_guide) == len(original)
        # markSilencers lowercases unchanged bases; an uppercase base marks a mutation.
        assert any(b.isupper() for b in silenced_guide)

    # A guide is either silenced or unsilenced, never both.
    assert set(silenced).isdisjoint(set(unsilenced))


def test_guide_silence_runs_on_unsilenceable_guides(fus3):
    _orf, orf_plus, start, end = fus3
    forward, reverse = L.findPamSites(orf_plus, start, end)
    guides = dict(forward)
    guides.update(reverse)
    _silenced, unsilenced = L.tryToPamSilence(orf_plus, start, end, guides)

    guide_silenced = L.guideSilence(orf_plus, start, unsilenced)
    # Result is keyed by a subset of the unsilenced guides.
    assert set(guide_silenced).issubset(set(unsilenced))
