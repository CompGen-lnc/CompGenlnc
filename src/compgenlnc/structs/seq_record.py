from __future__ import annotations
from dataclasses import dataclass
import copy

@dataclass(slots=True)
class SeqRecord:
    _seq: str

    def __init__(self, seq: str | SeqRecord, /):
        if isinstance(seq, SeqRecord):
            self = copy.deepcopy(seq)
            return

        self._seq = seq.strip()

    def __str__(self):
        return self._seq

    def __bool__(self):
        return bool(self._seq)

    def __len__(self):
        return len(self._seq)
    