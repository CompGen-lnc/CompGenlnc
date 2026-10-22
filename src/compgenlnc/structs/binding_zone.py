from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class BindingZone:
    _lnc_id: str
    _mir_id: str
    _hint: bool
    _energy: float
    _lnc_range: range
    _mir_range: range

    def __init__(
        self,
        lnc: str,
        mir: str,
        hint: bool,
        energy: float = 0,
        lnc_range: range = range(0),
        mir_range: range = range(0),
    ) -> None:
        self._lnc_id = lnc
        self._mir_id = mir
        self._hint = hint
        self._energy = energy
        self._lnc_range = lnc_range
        self._mir_range = mir_range

    @property
    def lncrna(self) -> str:
        return self._lnc_id

    @property
    def mirna(self) -> str:
        return self._mir_id

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def lnc_start(self) -> int:
        return self._lnc_range.start + 1

    @property
    def lnc_end(self) -> int:
        return self._lnc_range.stop

    @property
    def mir_start(self) -> int:
        return self._mir_range.start + 1

    @property
    def mir_end(self) -> int:
        return self._mir_range.stop

    def adjust_mirna_range(self, start: int) -> BindingZone:
        new_range = range(
            self._mir_range.start + start, self._mir_range.stop + start
        )
        return BindingZone(
            self._lnc_id,
            self._mir_id,
            self._hint,
            self._energy,
            self._lnc_range,
            new_range,
        )

    def __bool__(self) -> bool:
        return self._hint
