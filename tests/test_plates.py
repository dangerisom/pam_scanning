"""Unit tests for the vendored Microplate well-ordering.

The column-major fill order (A1, B1, ..., H1, A2, ...) was verified against
primer-order workbooks produced by the original code, so locking it here guards
byte-identical primer orders.
"""

import pytest

from pam_scanning.plates import Microplate


def test_96_well_is_column_major():
    plate = Microplate()
    plate.create(8, 12)
    wells = plate.plateArrayTransposed
    assert len(wells) == 96
    # First column fills A..H before moving to column 2.
    assert wells[:9] == ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2"]
    assert wells[-1] == "H12"


def test_384_well_dimensions():
    plate = Microplate()
    plate.create(16, 24)
    wells = plate.plateArrayTransposed
    assert len(wells) == 384
    assert wells[0] == "A1"
    assert wells[16] == "A2"  # 16 rows per column
    assert wells[-1] == "P24"


def test_create_returns_self_and_records_shape():
    plate = Microplate()
    assert plate.create(8, 12) is plate
    assert (plate.rows, plate.cols) == (8, 12)


def test_too_many_rows_raises():
    with pytest.raises(ValueError):
        Microplate().create(17, 24)
