"""Locate or install the external NCBI BLAST+ ``blastn`` executable.

PAM-scanning shells out to ``blastn`` for off-target evaluation, but BLAST+ is a
compiled C++ toolkit -- not a Python package -- so pip cannot provide it. Rather
than leave a user stuck at a bare ``[Errno 2] ... 'blastn'`` error, this module
can download the official NCBI BLAST+ binaries into an app-managed directory
(``~/.pam_scanning/blast``) and put just that directory on ``PATH`` for the run.

This is deliberately *non-intrusive*: nothing is installed into a conda
environment (so no base env gets weighed down), no conda is required at all, and
the user's shell profile is never edited. The binaries live under the app's own
directory and are re-activated on each run, so a lab member without conda -- or
anyone -- gets a working ``blastn`` with one click. It also drops in
``makeblastdb`` and the rest of the BLAST+ suite for building the local database.

The install is never automatic: callers confirm first (a dialog in the GUI, the
``--install-blast`` flag in the CLI).
"""

import glob
import hashlib
import os
import platform
import shutil
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request


_LATEST_URL = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
_BINARY = "blastn.exe" if os.name == "nt" else "blastn"


def app_blast_dir():
    """Directory where the app keeps its managed BLAST+ install."""
    return os.path.join(os.path.expanduser("~"), ".pam_scanning", "blast")


def which_blastn():
    """Return the path to ``blastn`` if it is on PATH, else ``None``."""
    return shutil.which("blastn")


def local_blastn():
    """Return the app-managed ``blastn`` if a prior download exists, else ``None``."""
    matches = glob.glob(os.path.join(app_blast_dir(), "**", _BINARY), recursive=True)
    return matches[0] if matches else None


def _add_to_path(directory):
    """Prepend *directory* to this process's PATH (idempotent)."""
    entries = os.environ.get("PATH", "").split(os.pathsep)
    if directory not in entries:
        os.environ["PATH"] = directory + os.pathsep + os.environ.get("PATH", "")


def activate_local_blast():
    """Put a previously downloaded app-managed BLAST+ on PATH; return its path or None."""
    found = local_blastn()
    if found:
        _add_to_path(os.path.dirname(found))
    return found


def ensure_available():
    """Return a usable ``blastn`` path (from PATH or a prior app install), or ``None``."""
    return which_blastn() or activate_local_blast()


# --- Platform selection ------------------------------------------------------

def _platform_tag():
    """NCBI filename fragment for this platform, or None if unsupported.

    Prefers a native build; on macOS falls back to the universal binary.
    """
    system, machine = platform.system(), platform.machine().lower()
    arm = machine in ("arm64", "aarch64")
    if system == "Darwin":
        return ["aarch64-macosx", "universal-macosx"] if arm else ["x64-macosx", "universal-macosx"]
    if system == "Linux":
        return ["aarch64-linux"] if arm else ["x64-linux"]
    if system == "Windows":
        return ["x64-win64"]
    return None


def _resolve_download(listing, tags):
    """Find the first ``ncbi-blast-<ver>+-<tag>.tar.gz`` present in a directory listing."""
    import re

    for tag in tags:
        match = re.search(r"ncbi-blast-[0-9.]+\+-%s\.tar\.gz" % re.escape(tag), listing)
        if match:
            return match.group(0)
    return None


# --- Download / install ------------------------------------------------------

def _http_get_text(url, timeout=60):
    request = urllib.request.Request(url, headers={"User-Agent": "pam_scanning/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")


def _download(url, dest, log, timeout=120):
    """Stream *url* to *dest*, logging coarse progress."""
    request = urllib.request.Request(url, headers={"User-Agent": "pam_scanning/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response, open(dest, "wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        got, last = 0, -10
        while True:
            chunk = response.read(1 << 16)
            if not chunk:
                break
            out.write(chunk)
            got += len(chunk)
            pct = int(got * 100 / total) if total else 0
            if total and pct >= last + 10:
                log("  downloaded %d%% (%d/%d MB)\n" % (pct, got >> 20, total >> 20))
                last = pct
    return dest


def _md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_md5(tarball, url, log):
    """Verify *tarball* against NCBI's ``.md5`` sidecar; warn (don't fail) if unavailable."""
    try:
        expected = _http_get_text(url + ".md5").split()[0]
    except (urllib.error.URLError, OSError, IndexError):
        log("  (no MD5 checksum available; skipping verification)\n")
        return
    if _md5(tarball) != expected:
        raise RuntimeError("Downloaded BLAST+ archive failed its MD5 check; aborting.")
    log("  MD5 verified.\n")


def _extract(tarball, into):
    """Extract *tarball* into *into*, using the safe 'data' filter where available."""
    with tarfile.open(tarball) as tar:
        try:
            tar.extractall(into, filter="data")   # Python >= 3.12
        except TypeError:                          # pragma: no cover - older Pythons
            tar.extractall(into)


def install_blast(log=print):
    """Download and activate NCBI BLAST+, returning the ``blastn`` path.

    *log* receives progress text as it streams. Raises :class:`RuntimeError` with
    an actionable message on any failure. If ``blastn`` is already available (on
    PATH or from a prior app install) it is returned without downloading.
    """
    found = ensure_available()
    if found:
        return found

    tags = _platform_tag()
    if tags is None:
        raise RuntimeError(
            "Automatic BLAST+ install is not supported on this platform (%s/%s). "
            "Install BLAST+ manually: conda install -c bioconda blast, or from "
            "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
            % (platform.system(), platform.machine()))

    log("Locating the latest NCBI BLAST+ build...\n")
    try:
        listing = _http_get_text(_LATEST_URL)
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError("Could not reach the NCBI download server (%s). Try again later, "
                           "or install BLAST+ manually." % getattr(exc, "reason", exc))
    filename = _resolve_download(listing, tags)
    if filename is None:
        raise RuntimeError("No NCBI BLAST+ build found for this platform (%s)." % ", ".join(tags))

    dest_dir = app_blast_dir()
    os.makedirs(dest_dir, exist_ok=True)
    url = _LATEST_URL + filename
    log("Downloading %s ...\n" % filename)
    fd, tmp = tempfile.mkstemp(suffix=".tar.gz", dir=dest_dir)
    os.close(fd)
    try:
        _download(url, tmp, log)
        _verify_md5(tmp, url, log)
        log("Extracting...\n")
        _extract(tmp, dest_dir)
    except (urllib.error.URLError, OSError, tarfile.TarError) as exc:
        raise RuntimeError("Failed to download/extract BLAST+: %s" % exc)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    path = activate_local_blast()
    if path is None:
        raise RuntimeError("BLAST+ was downloaded but 'blastn' could not be located afterward.")
    log("\nBLAST+ ready: %s\n" % path)
    return path


if __name__ == "__main__":   # convenience: `python -m pam_scanning.blast_setup`
    if ensure_available():
        print("blastn available:", which_blastn() or local_blastn())
    else:
        try:
            install_blast()
        except RuntimeError as exc:
            sys.exit(str(exc))
