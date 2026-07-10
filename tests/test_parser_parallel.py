"""Tests for the parallel BLAST-results parser plumbing.

The deep equivalence of per-block vs whole-file evaluation is verified against a
real BLAST run during development; these guard the splitting and dispatch logic
that make parallelization safe (no genome/BLAST needed).
"""

from pam_scanning import library


def test_split_guide_blocks_partitions_by_query():
    lines = [
        "BLASTN 2.17.0+\n", "\n",
        "Query= 1 AAAA guide query\n", "Length=4\n", "hit line\n", "\n",
        "Query= 2 CCCC guide query\n", "Length=4\n", "hit line\n", "\n",
    ]
    blocks = library._split_guide_blocks(lines)
    assert [b[0].split()[1] for b in blocks] == ["1", "2"]        # one block per guide
    # The header before the first Query= is not part of any block.
    assert all(not any("BLASTN" in line for line in b) for b in blocks)


def test_split_guide_blocks_empty_without_queries():
    assert library._split_guide_blocks(["no queries here\n", "\n"]) == []


def test_evaluate_blocks_serial_returns_one_result_per_block():
    # Two no-hit guides: exercises the serial dispatch + result shape without any
    # genome scanning. A no-hit guide contributes nothing, mirroring the loop.
    b1 = ["Query= 1 " + "A" * 23 + " guide query\n", "***** No hits found *****\n", "\n"]
    b2 = ["Query= 2 " + "C" * 23 + " guide query\n", "***** No hits found *****\n", "\n"]
    results = library._evaluate_blocks([b1, b2], {}, "unused-path", 15, 5)
    assert len(results) == 2
    for safe, unsafe, pam_incl, super_cons in results:
        assert (safe, unsafe, pam_incl, super_cons) == ({}, {}, {}, {})


def test_parse_blast_lines_whole_file_equals_split(monkeypatch):
    # Whatever _parse_blast_lines yields over the whole file must equal the merge
    # of per-block results -- the invariant the parallel path relies on.
    lines = [
        "Query= 1 " + "A" * 23 + " guide query\n", "***** No hits found *****\n", "\n",
        "Query= 2 " + "C" * 23 + " guide query\n", "***** No hits found *****\n", "\n",
    ]
    whole = library._parse_blast_lines(lines, {}, 15, 5)
    merged = ({}, {}, {}, {})
    for block in library._split_guide_blocks(lines):
        for acc, part in zip(merged, library._parse_blast_lines(block, {}, 15, 5)):
            acc.update(part)
    assert merged == whole
