from __future__ import annotations

import copy
from dataclasses import dataclass


@dataclass(slots=True)
class StructureRecord:
    _ss: str

    def __init__(self, ss: str | StructureRecord, /):
        if isinstance(ss, StructureRecord):
            self._ss = copy.deepcopy(ss._ss)
            return

        self._ss = ss.strip()

    def __str__(self):
        return self._ss

    def __bool__(self):
        return bool(self._ss)

    def __len__(self):
        return len(self._ss)

    def __getitem__(self, key):
        return self._ss[key]

    def __iter__(self):
        return iter(self._ss)

    def __eq__(self, other):
        if isinstance(other, str):
            return self._ss == other
        elif not isinstance(other, StructureRecord):
            return NotImplemented
        return self._ss == other._ss
