"""Test the bundled default yeast genome resolver (decompress + cache)."""

import os

from pam_scanning import chimeras


def test_default_genome_decompresses_and_caches(tmp_path, monkeypatch):
    # Redirect the cache to a temp home so the test never touches the real one.
    monkeypatch.setattr(os.path, "expanduser",
                        lambda p: str(tmp_path) if p == "~" else os.path.realpath(p))

    path = chimeras.default_genome_path()
    assert os.path.isfile(path)
    assert path.endswith(chimeras.DEFAULT_GENOME_NAME)
    assert os.path.getsize(path) > 10_000_000          # ~12 MB expanded
    with open(path) as fh:
        assert fh.readline().startswith(">")           # FASTA header

    # A second call returns the same cached file without re-expanding.
    mtime = os.path.getmtime(path)
    assert chimeras.default_genome_path() == path
    assert os.path.getmtime(path) == mtime
