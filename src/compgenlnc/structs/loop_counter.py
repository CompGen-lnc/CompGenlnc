from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np


loop_tuple = np.dtype(
    [
        ("kind", "U10"),
        ("low", np.uint16),
        ("high", np.uint16),
        ("energy", np.int16),
    ]
)


@dataclass(slots=True)
class LoopCounter:
    _energy: float
    _loops: np.typing.NDArray[np.void]
    _hairpin_energy: int
    _interior_energy: int
    _multi_energy: int

    def __init__(
        self,
        loops: Iterable[tuple[str, int, int, int]],
        energy: float | None = None,
    ) -> None:
        loops_arr = np.array(loops, loop_tuple)
        self._loops = loops_arr[loops_arr["kind"] != "External"]
        if energy is not None:
            self._energy = energy
        else:
            self._energy = self._loops["energy"].sum()
        self._hairpin_energy = self.hairpin_loops["energy"].sum()
        self._interior_energy = self.interior_loops["energy"].sum()
        self._multi_energy = self.multi_loops["energy"].sum()

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def loops(self) -> np.typing.NDArray[np.void]:
        return self._loops

    @property
    def hairpin_loops(self) -> np.typing.NDArray[np.void]:
        mask = self._loops["kind"] == "Hairpin"
        return self._loops[mask]

    @property
    def hairpin_energy(self) -> int:
        return self._hairpin_energy

    @property
    def interior_loops(self) -> np.typing.NDArray[np.void]:
        mask = self._loops["kind"] == "Interior"
        return self._loops[mask]

    @property
    def interior_energy(self) -> int:
        return self._interior_energy

    @property
    def multi_loops(self) -> np.typing.NDArray[np.void]:
        mask = self._loops["kind"] == "Multi"
        return self._loops[mask]

    @property
    def multi_energy(self) -> int:
        return self._multi_energy
