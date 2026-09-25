from typing import Annotated

import numpy as np


dtype_loop_tuple = np.dtype(
    [
        ("kind", "U10"),
        ("low", np.uint16),
        ("high", np.uint16),
        ("energy", np.int16),
    ]
)
dtype_interaction_tuple = np.dtype(
    [
        ("lncRNA", np.dtypes.StringDType),
        ("miRNA", np.dtypes.StringDType),
        ("positive", np.bool),
    ]
)

dtype_interaction_pair = np.dtype(
    [
        ("lncRNA", np.dtypes.StringDType),
        ("miRNA", np.dtypes.StringDType),
    ]
)

type LoopTuple = Annotated[np.void, dtype_loop_tuple]
type InteractionTuple = Annotated[np.void, dtype_interaction_tuple]
type InteractionPair = (
    Annotated[np.void, dtype_interaction_pair] | InteractionTuple
)
