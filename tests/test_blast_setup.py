"""Tests for BLAST+ discovery/activation logic (no network, no download).

Only pure logic is exercised: platform selection, parsing an NCBI directory
listing, locating an app-managed install, and activating it onto PATH. The
actual download is never performed.
"""

import os

import pytest

from pam_scanning import blast_setup as bs


# --- Platform selection ----------------------------------------------------

def test_platform_tag_is_known_here():
    # On any of the supported OSes this must resolve to at least one build tag.
    tags = bs._platform_tag()
    assert tags is None or (isinstance(tags, list) and tags)


def test_platform_tag_macos_arm(monkeypatch):
    monkeypatch.setattr(bs.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(bs.platform, "machine", lambda: "arm64")
    assert bs._platform_tag()[0] == "aarch64-macosx"


def test_platform_tag_linux_x64(monkeypatch):
    monkeypatch.setattr(bs.platform, "system", lambda: "Linux")
    monkeypatch.setattr(bs.platform, "machine", lambda: "x86_64")
    assert bs._platform_tag() == ["x64-linux"]


# --- Parsing the NCBI directory listing ------------------------------------

_LISTING = (
    'ncbi-blast-2.17.0+-aarch64-linux.tar.gz  ncbi-blast-2.17.0+-aarch64-macosx.tar.gz '
    'ncbi-blast-2.17.0+-universal-macosx.tar.gz  ncbi-blast-2.17.0+-x64-linux.tar.gz'
)


def test_resolve_download_picks_first_matching_tag():
    got = bs._resolve_download(_LISTING, ["aarch64-macosx", "universal-macosx"])
    assert got == "ncbi-blast-2.17.0+-aarch64-macosx.tar.gz"


def test_resolve_download_falls_back_to_universal():
    got = bs._resolve_download(_LISTING, ["x64-macosx", "universal-macosx"])
    assert got == "ncbi-blast-2.17.0+-universal-macosx.tar.gz"


def test_resolve_download_none_when_unavailable():
    assert bs._resolve_download(_LISTING, ["x64-win64"]) is None


# --- App-managed install location + activation -----------------------------

def test_local_blastn_none_when_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(bs, "app_blast_dir", lambda: str(tmp_path))
    assert bs.local_blastn() is None


def test_activate_local_blast_adds_to_path(tmp_path, monkeypatch):
    bindir = tmp_path / "ncbi-blast-2.17.0+" / "bin"
    bindir.mkdir(parents=True)
    (bindir / bs._BINARY).write_text("#!/bin/sh\n")
    monkeypatch.setattr(bs, "app_blast_dir", lambda: str(tmp_path))
    monkeypatch.setenv("PATH", "/usr/bin")

    path = bs.activate_local_blast()
    assert path == str(bindir / bs._BINARY)
    assert str(bindir) in os.environ["PATH"].split(os.pathsep)


def test_ensure_available_prefers_path(monkeypatch):
    monkeypatch.setattr(bs, "which_blastn", lambda: "/usr/local/bin/blastn")
    # local_blastn must not be consulted when PATH already has it.
    monkeypatch.setattr(bs, "local_blastn", lambda: (_ for _ in ()).throw(AssertionError))
    assert bs.ensure_available() == "/usr/local/bin/blastn"


def test_install_blast_returns_existing_without_downloading(monkeypatch):
    monkeypatch.setattr(bs, "which_blastn", lambda: "/usr/local/bin/blastn")
    # _platform_tag / network must never be reached.
    monkeypatch.setattr(bs, "_platform_tag", lambda: (_ for _ in ()).throw(AssertionError))
    assert bs.install_blast(log=lambda *_: None) == "/usr/local/bin/blastn"


# --- Building a BLAST database from a genome --------------------------------

def test_blast_db_exists(tmp_path):
    prefix = str(tmp_path / "yeast")
    assert not bs.blast_db_exists(prefix)
    open(prefix + ".nin", "w").close()
    assert bs.blast_db_exists(prefix)


def test_genome_key_is_stable_and_size_sensitive(tmp_path):
    g = tmp_path / "genome.fa"
    g.write_text(">chr1\nACGT\n")
    key = bs._genome_key(str(g))
    assert key == bs._genome_key(str(g))          # same file -> same key
    g.write_text(">chr1\nACGTACGT\n")             # content/size changed
    assert bs._genome_key(str(g)) != key          # -> rebuild


def test_ensure_blast_db_reuses_cached_without_rebuilding(tmp_path, monkeypatch):
    g = tmp_path / "genome.fa"
    g.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(bs, "app_blastdb_dir", lambda: str(tmp_path / "blastdb"))
    monkeypatch.setattr(bs, "blast_db_exists", lambda prefix: True)   # pretend it's built
    monkeypatch.setattr(bs, "ensure_available", lambda: (_ for _ in ()).throw(AssertionError))
    prefix = bs.ensure_blast_db(str(g), log=lambda *_: None)
    assert prefix.endswith("genome")             # <cache>/<genome-basename>


def test_ensure_blast_db_errors_without_makeblastdb(tmp_path, monkeypatch):
    g = tmp_path / "genome.fa"
    g.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(bs, "app_blastdb_dir", lambda: str(tmp_path / "blastdb"))
    monkeypatch.setattr(bs, "blast_db_exists", lambda prefix: False)
    monkeypatch.setattr(bs, "ensure_available", lambda: None)
    monkeypatch.setattr(bs.shutil, "which", lambda name: None)        # no makeblastdb
    with pytest.raises(RuntimeError, match="makeblastdb"):
        bs.ensure_blast_db(str(g), log=lambda *_: None)
