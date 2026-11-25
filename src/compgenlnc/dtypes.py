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
"""Data type for database indexes.

It has the following fields:

- `LNCRNA_FIELD` : *str*
- `PRECURSOR_FIELD` : *str*
- `MATURE_FIELD` : *str*
- `INTERACTION_FIELD` : *str*
"""
dtype_interaction_pair = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
    ]
)
"""Data type for interaction RNA ID pairs.

It has the following fields:

- `LNCRNA_FIELD` : *str*
- `MIRNA_FIELD` : *str*
"""
dtype_interaction_tuple = np.dtype(
    [
        (fields.LNCRNA_FIELD, np.dtypes.StringDType),
        (fields.MIRNA_FIELD, np.dtypes.StringDType),
        (fields.IS_POSITIVE_FIELD, np.bool),
    ]
)
"""Data type for marked interaction pairs.

It has the following fields:

- `LNCRNA_FIELD` : *str*
- `MIRNA_FIELD` : *str*
- `IS_POSITIVE_FIELD` : *bool*
"""
dtype_loop_tuple = np.dtype(
    [
        (fields.KIND_FIELD, np.dtypes.StringDType),
        (fields.LOW_POSITION, np.uint16),
        (fields.HIGH_POSITION, np.uint16),
        (fields.ENERGY_FIELD, np.int16),
    ]
)
"""Data type for loops tuples.

It has the following fields:

- `KIND_FIELD` : *str*
- `LOW_POSITION` : *int*
- `HIGH_POSITION` : *int*
- `ENERGY_FIELD` : *int*
"""
