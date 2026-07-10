"""Tests for the GUI's last-directory persistence (no display required).

Importing :mod:`pam_scanning.gui` only needs the tkinter module to import, not a
running display, so these run headless. Only the pure state helpers are exercised.
"""

import os

from pam_scanning import gui


def test_last_dir_round_trips(tmp_path, monkeypatch):
    state = tmp_path / "nested" / "gui_state.json"
    monkeypatch.setattr(gui, "_STATE_PATH", state)
    gui._save_last_dir(str(tmp_path))
    assert gui._load_last_dir() == str(tmp_path)
    assert state.is_file()  # parent directory was created


def test_missing_state_falls_back_to_cwd(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "_STATE_PATH", tmp_path / "absent.json")
    assert gui._load_last_dir() == os.getcwd()


def test_stale_directory_is_ignored(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "_STATE_PATH", tmp_path / "gui_state.json")
    gui._save_last_dir(str(tmp_path / "was_deleted"))  # never existed
    assert gui._load_last_dir() == os.getcwd()


def test_corrupt_state_is_tolerated(tmp_path, monkeypatch):
    state = tmp_path / "gui_state.json"
    state.write_text("{ not valid json")
    monkeypatch.setattr(gui, "_STATE_PATH", state)
    assert gui._load_last_dir() == os.getcwd()  # no exception
