from typing import Iterable

import numpy as np

from compgenlnc.structs import SeqSSRecord
from compgenlnc.typing import SSLike


def recode_2d_structure(ss: SSLike) -> np.typing.NDArray[np.uint8]:
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
    ss_list: Iterable[SSLike],
) -> Iterable[np.typing.NDArray[np.uint8]]:
    return (recode_2d_structure(ss) for ss in ss_list)
