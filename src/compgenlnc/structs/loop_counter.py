from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from compgenlnc.consts.fields import ENERGY_FIELD, KIND_FIELD
from compgenlnc.dtypes import dtype_loop_tuple
from compgenlnc.typing.numpy_dtypes import LoopTuple


@dataclass(slots=True)
class LoopCounter:
    """Container for the loops of a 2D structure.

    Attributes
    ----------
    rna : str
        The ID of the RNA.
    energy : float
        Total MFE from the loops.
    loops : numpy.ndarray of LoopTuple
        A 1D array with the loops.
    hairpin_loops : numpy.ndarray of LoopTuple
        A 1D array with the hairpin loops.
    hairpin_energy : float
        Total MFE from the hairpin loops.
    stack_loops : numpy.ndarray of LoopTuple
        A 1D array with the stack loops.
    stack_energy : float
        Total MFE from the stack loops.
    multi_loops : numpy.ndarray of LoopTuple
        A 1D array with the multi loops.
    multi_energy : float
        Total MFE from the multi loops.
    loops_vectorized : numpy.ndarray of float32
        A 1D array with the energy and count of each loop kind.
    """

    _rna: str
    _energy: float
    _loops: NDArray[np.void]
    _hairpin_energy: float
    _stack_energy: float
    _multi_energy: float

    def __init__(
        self,
        id_: str,
        loops: Sequence[LoopTuple],
        energy: float | None = None,
    ) -> None:
        """Initialize a LoopCounter object.

        Parameters
        ----------
        id_: str
            The ID of the RNA.
        loops: sequence of LoopTuple
            The loops of the RNA.
        energy : float, optional
            The energy of the RNA. If not provided, it will be
            calculated from the loops. Default: `None`.
        """
        self._rna = id_
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
    def rna(self) -> str:
        """The ID of the RNA."""
        return self._rna

    @property
    def energy(self) -> float:
        """Total MFE from the loops."""
        return self._energy

    @property
    def loops(self) -> NDArray[LoopTuple]:
        """Array with the loops."""
        return self._loops

    @property
    def hairpin_loops(self) -> NDArray[LoopTuple]:
        """Array with the hairpin loops."""
        mask = self._loops[KIND_FIELD] == "Hairpin"
        return self._loops[mask]

    @property
    def hairpin_energy(self) -> int:
        """Total MFE from the hairpin loops."""
        return self._hairpin_energy

    @property
    def stack_loops(self) -> NDArray[LoopTuple]:
        """Array with the stack loops."""
        mask = self._loops[KIND_FIELD] == "Interior"
        return self._loops[mask]

    @property
    def stack_energy(self) -> int:
        """Total MFE from the stack loops."""
        return self._stack_energy

    @property
    def multi_loops(self) -> NDArray[LoopTuple]:
        """Array with the multi loops."""
        mask = self._loops[KIND_FIELD] == "Multi"
        return self._loops[mask]

    @property
    def multi_energy(self) -> int:
        """Total MFE from the multi loops."""
        return self._multi_energy

    @property
    def loops_vectorized(self) -> NDArray[np.float32]:
        """Array with the energy and count of each loop kind."""
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

    def __eq__(self, value: LoopCounter) -> bool:
        """Compare two LoopCounter objects.

        Parameters
        ----------
        value : LoopCounter
            The LoopCounter object to compare with.

        Returns
        -------
        bool
            `True` if the energy and loops are the same, `False`
            otherwise.
        """
        same_energy = self.energy == value.energy
        same_loops = np.array_equal(self.loops, value.loops)
        return same_energy and same_loops

    def __iter__(self) -> Iterator[LoopTuple]:
        """Iterate over the loops.

        Returns
        -------
        Iterator of LoopTuple
            An iterator over the loops.
        """
        return iter(self.loops)
