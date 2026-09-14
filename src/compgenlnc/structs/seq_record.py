from __future__ import annotations

import copy
from dataclasses import dataclass


@dataclass(slots=True)
class SeqRecord:
    _seq: str

    def __init__(self, seq: str | SeqRecord, /):
        if isinstance(seq, SeqRecord):
            self._seq = copy.deepcopy(seq._seq)
            return

        self._seq = seq.strip()

    def __str__(self):
        return self._seq

    def __bool__(self):
        return bool(self._seq)

    def __len__(self):
        return len(self._seq)

    def __getitem__(self, key):
        return self._seq[key]

    def __iter__(self):
        return iter(self._seq)

    def __eq__(self, other):
        if isinstance(other, str):
            return self._seq == other
        elif not isinstance(other, SeqRecord):
            return NotImplemented
        return self._seq == other._seq

    def convert_to_T(self):
        return SeqRecord(self._seq.replace("U", "T"))

    def convert_to_U(self):
        return SeqRecord(self._seq.replace("T", "U"))
