from collections.abc import Iterable

import numpy as np
from numpy.typing import NDArray

from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SSLike


def recode_2d_structure(ss: SSLike, /) -> NDArray[np.uint8]:
    """Recode a 2D structure into a numerical vector.

    Parameters
    ----------
    ss : SSLike
        The 2D structure string on dot-bracket notation.

    Returns
    -------
    numpy.ndarray of uint8
        A 1D array with the recoded structure.
    """
    if isinstance(ss, SeqSSRecord):
        ss = ss.ss

    recode_map = {".": 5, "(": 6, ")": 7}
    recoded = np.fromiter(
        (recode_map[char] for char in ss if char in recode_map.keys()),
        np.uint8,
        len(ss),
    )

    return recoded


def recode_2d_structure_list(
    structures: Iterable[SSLike], /
) -> Iterable[NDArray[np.uint8]]:
    """Recode a list of structures into numerical vectors.

    Parameters
    ----------
    structures : iterable object of SSLike
        Group of 2D structure strings on dot-bracket notation.

    Yields
    ------
    numpy.ndarray of uint8
        A 1D array with the recoded structure.
    """
    return (recode_2d_structure(ss) for ss in structures)
