from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class BindingZone:
    _hint: bool
    _energy: float
    _lnc_range: range
    _mir_range: range

    def __init__(self, hint: bool, energy: float = 0, lnc_range: range = range(0), mir_range: range = range(0)) -> None:
        self._hint = hint
        self._energy = energy
        self._lnc_range = lnc_range
        self._mir_range = mir_range

    @property
    def hint(self) -> bool:
        return self._hint

    @property
    def energy(self) -> bool:
        return self._energy

    @property
    def lnc_start(self) -> bool:
        return self._lnc_range.start + 1

    @property
    def lnc_end(self) -> bool:
        return self._lnc_range.stop

    @property
    def mir_start(self) -> bool:
        return self._mir_range.start + 1

    @property
    def mir_end(self) -> bool:
        return self._mir_range.stop
