from typing import Annotated

import numpy as np

import compgenlnc.consts.fields as fields


dtype_loop_tuple = np.dtype(
    [
        (fields.KIND_FIELD, "U10"),
        (fields.LOW_POSITION, np.uint16),
        (fields.HIGH_POSITION, np.uint16),
        (fields.ENERGY_FIELD, np.int16),
    ]
)
dtype_interaction_tuple = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
        (fields.IS_POSITIVE_FIELD, np.bool),
    ]
)
dtype_interaction_pair = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
    ]
)
dtype_index = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.PRECURSOR_FIELD, np.dtypes.StringDType),
        (fields.MATURE_FIELD, np.dtypes.StringDType),
    ]
)

type LoopTuple = Annotated[np.void, dtype_loop_tuple]
type InteractionTuple = Annotated[np.void, dtype_interaction_tuple]
type InteractionPair = (
    Annotated[np.void, dtype_interaction_pair] | InteractionTuple
)
type Index = Annotated[np.void, dtype_index]
