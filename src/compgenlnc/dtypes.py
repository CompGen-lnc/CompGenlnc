import numpy as np

import compgenlnc.consts.fields as fields

dtype_index = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.PRECURSOR_FIELD, np.dtypes.StringDType),
        (fields.MATURE_FIELD, np.dtypes.StringDType),
        (fields.INTERACTION_FIELD, np.dtypes.StringDType),
    ]
)
dtype_interaction_pair = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
    ]
)
dtype_interaction_tuple = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
        (fields.IS_POSITIVE_FIELD, np.bool),
    ]
)
dtype_loop_tuple = np.dtype(
    [
        (fields.KIND_FIELD, np.dtypes.StringDType),
        (fields.LOW_POSITION, np.uint16),
        (fields.HIGH_POSITION, np.uint16),
        (fields.ENERGY_FIELD, np.int16),
    ]
)
