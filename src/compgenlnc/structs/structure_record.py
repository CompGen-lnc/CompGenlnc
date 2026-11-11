from __future__ import annotations

import copy
from collections.abc import Iterator
from dataclasses import dataclass


@dataclass(slots=True)
class StructureRecord:
    """Record for the RNA 2D structure."""

    _ss: str

    def __init__(self, ss: str | StructureRecord, /) -> None:
        """Initialize the StructureRecord with a sequence.

        Parameters
        ----------
        ss : str or StructureRecord
            The 2D structure to initialize the StructureRecord with. If
            a StructureRecord is provided, the 2D structure will be
            copied.
        """
        if isinstance(ss, StructureRecord):
            self._ss = copy.deepcopy(ss._ss)
            return

        self._ss = ss.strip()

    def __str__(self) -> str:
        """Get the 2D structure as a string.

        Returns
        -------
        str
            The 2D structure as a string.
        """
        return self._ss

    def __bool__(self) -> bool:
        """Whether the 2D structure is not empty.

        Returns
        -------
        bool
            `True` if the 2D structure is not empty, `False` otherwise.
        """
        return bool(self._ss)

    def __len__(self) -> int:
        """Get the length of the 2D structure.

        Returns
        -------
        int
            The length of the 2D structure.
        """
        return len(self._ss)

    def __getitem__(self, key) -> str:
        """Get a substring from the 2D structure by index.

        Parameters
        ----------
        key : int or slice of int
           The index or slice of the 2D structure to retrieve.

        Returns
        -------
        str
            The substring from the 2D structure.
        """
        return self._ss[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate over the 2D structure.

        Returns
        -------
        Iterator of str
            An iterator over the 2D structure.
        """
        return iter(self._ss)

    def __eq__(self, other) -> bool:
        """Whether the 2D structure is equal to another one.

        Parameters
        ----------
        other : str or SeqRecord
            The 2D structure to compare with.

        Returns
        -------
        bool
            `True` if the 2D structures are equal, `False` otherwise.
        """
        if isinstance(other, str):
            return self._ss == other
        elif not isinstance(other, StructureRecord):
            return NotImplemented
        return self._ss == other._ss
