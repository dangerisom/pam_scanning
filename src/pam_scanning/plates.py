"""Minimal microplate helper used for laying out primer orders.

This is a self-contained replacement for the previously external (and missing)
``biomek.platesBiomek.Microplate`` dependency. Only the surface actually used by
:func:`pam_scanning.library.createPrimerOrder` is implemented:

* :meth:`Microplate.create` builds the plate of a given size.
* :attr:`Microplate.plateArrayTransposed` is the list of well labels in
  **column-major** order (A1, B1, ..., H1, A2, ...), i.e. filled down each column
  before moving to the next. This ordering was verified against primer-order
  workbooks produced by the original code, so generated orders are byte-identical.
"""

# Standard microtiter-plate row letters (A-P covers up to 384-well plates).
_ROW_LETTERS = "ABCDEFGHIJKLMNOP"


class Microplate:
    """A rectangular microtiter plate addressed by well labels such as ``"A1"``."""

    def __init__(self):
        self.rows = 0
        self.cols = 0
        #: Well labels in column-major fill order (A1, B1, ..., A2, B2, ...).
        self.plateArrayTransposed = []

    def create(self, rows, cols):
        """Create a ``rows`` x ``cols`` plate (e.g. 8x12 = 96, 16x24 = 384)."""
        if rows > len(_ROW_LETTERS):
            raise ValueError(
                "Plate has %d rows but only %d row letters are defined"
                % (rows, len(_ROW_LETTERS))
            )
        letters = _ROW_LETTERS[:rows]
        self.rows = rows
        self.cols = cols
        self.plateArrayTransposed = [
            "%s%d" % (letters[r], c + 1)
            for c in range(cols)
            for r in range(rows)
        ]
        return self
