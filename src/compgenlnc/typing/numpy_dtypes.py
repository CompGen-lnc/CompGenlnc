from typing import Annotated

import numpy as np

from compgenlnc.dtypes import (
    dtype_index,
    dtype_interaction_pair,
    dtype_interaction_tuple,
    dtype_loop_tuple,
)

type Index = Annotated[np.void, dtype_index]
type InteractionTuple = Annotated[np.void, dtype_interaction_tuple]
type InteractionPair = (
    Annotated[np.void, dtype_interaction_pair] | InteractionTuple
)
type LoopTuple = Annotated[np.void, dtype_loop_tuple]
