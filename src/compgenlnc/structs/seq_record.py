from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import SupportsIndex


@dataclass(slots=True)
class SeqRecord:
    _seq: str

    def __init__(self, seq: str | SeqRecord, /) -> SeqRecord:
        if isinstance(seq, SeqRecord):
            self._seq = copy.deepcopy(seq._seq)
            return

        self._seq = seq.strip()

    def __str__(self) -> str:
        return self._seq

    def __bool__(self) -> bool:
        return bool(self._seq)

    def __len__(self) -> int:
        return len(self._seq)

    def __getitem__(
        self, key: SupportsIndex | slice[SupportsIndex | None]
    ) -> str:
        return self._seq[key]

    def __iter__(self) -> str:
        return iter(self._seq)

    def __contains__(self, item: str | SeqRecord) -> bool:
        if hasattr(item, "seq"):
            item = item.seq
        return str(item).upper() in self._seq.upper()

    def __eq__(self, other: str | SeqRecord) -> bool:
        if isinstance(other, str):
            return self._seq == other
        elif not isinstance(other, SeqRecord):
            return NotImplemented
        return self._seq == other._seq

    def convert_to_T(self) -> SeqRecord:
        return SeqRecord(self._seq.replace("U", "T"))

    def convert_to_U(self) -> SeqRecord:
        return SeqRecord(self._seq.replace("T", "U"))
