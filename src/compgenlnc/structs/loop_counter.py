from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

import numpy as np

from compgenlnc.consts.fields import ENERGY_FIELD, KIND_FIELD
from compgenlnc.typing.numpy_dtypes import dtype_loop_tuple, LoopTuple


@dataclass(slots=True)
class LoopCounter:
    _energy: float
    _loops: np.typing.NDArray[np.void]
    _hairpin_energy: float
    _stack_energy: float
    _multi_energy: float

    def __init__(
        self,
        loops: Iterable[tuple[str, int, int, int]],
        energy: float | None = None,
    ) -> None:
        loops_arr = np.array(loops, dtype_loop_tuple)
        self._loops = loops_arr
        if energy is not None:
            self._energy = energy
        else:
            self._energy = self._loops[ENERGY_FIELD].sum() / 100
        self._hairpin_energy = self.hairpin_loops[ENERGY_FIELD].sum() / 100
        self._stack_energy = self.stack_loops[ENERGY_FIELD].sum() / 100
        self._multi_energy = self.multi_loops[ENERGY_FIELD].sum() / 100

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def loops(self) -> np.typing.NDArray[LoopTuple]:
        return self._loops

    @property
    def hairpin_loops(self) -> np.typing.NDArray[LoopTuple]:
        mask = self._loops[KIND_FIELD] == "Hairpin"
        return self._loops[mask]

    @property
    def hairpin_energy(self) -> int:
        return self._hairpin_energy

    @property
    def stack_loops(self) -> np.typing.NDArray[LoopTuple]:
        mask = self._loops[KIND_FIELD] == "Interior"
        return self._loops[mask]

    @property
    def stack_energy(self) -> int:
        return self._stack_energy

    @property
    def multi_loops(self) -> np.typing.NDArray[LoopTuple]:
        mask = self._loops[KIND_FIELD] == "Multi"
        return self._loops[mask]

    @property
    def multi_energy(self) -> int:
        return self._multi_energy

    @property
    def loops_vectorized(self) -> np.typing.NDArray[np.float32]:
        return np.array(
            [
                len(self.hairpin_loops),
                len(self.stack_loops),
                len(self.multi_loops),
                self.energy,
                self.hairpin_energy,
                self.stack_energy,
                self.multi_energy,
            ],
            np.float32,
        )

    def __eq__(self, value: LoopCounter):
        same_energy = self.energy == value.energy
        same_loops = np.array_equal(self.loops, value.loops)
        return same_energy and same_loops

    def __iter__(self) -> Iterator[LoopTuple]:
        return iter(self.loops)
