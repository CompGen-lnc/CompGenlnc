from typing import Annotated

import numpy as np


loop_tuple = np.dtype(
    [
        ("kind", "U10"),
        ("low", np.uint16),
        ("high", np.uint16),
        ("energy", np.int16),
    ]
)

type LoopTuple = Annotated[np.void, loop_tuple]