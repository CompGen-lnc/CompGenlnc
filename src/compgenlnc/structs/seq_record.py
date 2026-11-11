from __future__ import annotations

import copy
from collections.abc import Iterator
from dataclasses import dataclass
from typing import SupportsIndex


@dataclass(slots=True)
class SeqRecord:
    """Record for the RNA sequence."""

    _seq: str

    def __init__(self, seq: str | SeqRecord, /) -> None:
        """Initialize the SeqRecord with a sequence.

        Parameters
        ----------
        seq : str or SeqRecord
            The sequence to initialize the SeqRecord with. If a
            SeqRecord is provided, the sequence will be copied.
        """
        if isinstance(seq, SeqRecord):
            self._seq = copy.deepcopy(seq._seq)
            return

        self._seq = seq.strip()

    def __str__(self) -> str:
        """Get the sequence as a string.

        Returns
        -------
        str
            The sequence as a string.
        """
        return self._seq

    def __bool__(self) -> bool:
        """Verify whether the sequence is not empty.

        Returns
        -------
        bool
            `True` if the sequence is not empty, `False` otherwise.
        """
        return bool(self._seq)

    def __len__(self) -> int:
        """Get the length of the sequence.

        Returns
        -------
        int
            The length of the sequence.
        """
        return len(self._seq)

    def __getitem__(
        self, key: SupportsIndex | slice[SupportsIndex | None]
    ) -> str:
        """Get a substring from the sequence by index.

        Parameters
        ----------
        key : int or slice of int
           The index or slice of the sequence to retrieve.

        Returns
        -------
        str
            The substring from the sequence.
        """
        return self._seq[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate over the sequence nucleotides.

        Returns
        -------
        Iterator of str
            An iterator over the sequence nucleotides.
        """
        return iter(self._seq)

    def __contains__(self, item: str | SeqRecord) -> bool:
        """Verify whether the sequence contains another one.

        Parameters
        ----------
        item : str or SeqRecord
            The sequence to check for containment.

        Returns
        -------
        bool
            `True` if the sequence contains the other one, `False`
            otherwise.
        """
        if hasattr(item, "__str__"):
            item = str(item)
        return str(item).upper() in self._seq.upper()

    def __eq__(self, other: str | SeqRecord) -> bool:
        """Whether the sequence is equal to another one.

        Parameters
        ----------
        other : str or SeqRecord
            The sequence to compare with.

        Returns
        -------
        bool
            `True` if the sequences are equal, `False` otherwise.
        """
        if isinstance(other, str):
            return self._seq == other
        elif not isinstance(other, SeqRecord):
            return NotImplemented
        return self._seq == other._seq

    def convert_to_T(self) -> SeqRecord:
        """Convert the sequence to ATCG alphabet.

        Returns
        -------
        SeqRecord
            The sequence converted to ATCG alphabet.
        """
        return SeqRecord(self._seq.replace("U", "T"))

    def convert_to_U(self) -> SeqRecord:
        """Convert the sequence to AUCG alphabet.

        Returns
        -------
        SeqRecord
            The sequence converted to AUCG alphabet.
        """
        return SeqRecord(self._seq.replace("T", "U"))
