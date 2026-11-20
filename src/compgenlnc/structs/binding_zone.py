from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class BindingZone:
    """Description for a binding zone of an interaction.

    Attributes
    ----------
    lncrna : str
        The ID of the interaction lncRNA.
    mirna : str
        The ID of the interaction miRNA.
    energy : float
        The energy of the interaction.
    lnc_start : int
        The start position of the interaction lncRNA.
    lnc_end : int
        The end position of the interaction lncRNA.
    mir_start : int
        The start position of the interaction miRNA.
    mir_end : int
        The end position of the interaction miRNA.
    """

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
        """Initialize a BindingZone object.

        Parameters
        ----------
        lnc : str
            lncRNA ID.
        mir : str
            miRNA ID.
        hint : bool
            Whether a binding zone was found for the interaction.
        energy : float, optional
            The energy of the binding zone. Default: `0`.
        lnc_range : range, optional
            The range of the lncRNA that is in the binding zone.
            Default: `range(0)`
        mir_range : range, optional
            The range of the miRNA that is in the binding zone.
            Default: `range(0)`
        """
        self._lnc_id = lnc
        self._mir_id = mir
        self._hint = hint
        self._energy = energy
        self._lnc_range = lnc_range
        self._mir_range = mir_range

    @property
    def lncrna(self) -> str:
        """The ID of the interaction lncRNA."""
        return self._lnc_id

    @property
    def mirna(self) -> str:
        """The ID of the interaction miRNA."""
        return self._mir_id

    @property
    def energy(self) -> float:
        """The energy of the interaction."""
        return self._energy

    @property
    def lnc_start(self) -> int:
        """The start position of the interaction lncRNA."""
        return self._lnc_range.start + 1

    @property
    def lnc_end(self) -> int:
        """The end position of the interaction lncRNA."""
        return self._lnc_range.stop

    @property
    def mir_start(self) -> int:
        """The start position of the interaction miRNA."""
        return self._mir_range.start + 1

    @property
    def mir_end(self) -> int:
        """The end position of the interaction miRNA."""
        return self._mir_range.stop

    def fit_mirna_range(self, start: int, /) -> BindingZone:
        """Fit the miRNA range for the precursor.

        Parameters
        ----------
        start : int
            The start index of the miRNA in the precursor.

        Returns
        -------
        BindingZone
            A new BindingZone object with the updated miRNA range.
        """
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
        """Verify whether a binding zone was found for the interaction.

        Returns
        -------
        bool
            `True` if the binding zone exists, `False` otherwise.
        """
        return self._hint
