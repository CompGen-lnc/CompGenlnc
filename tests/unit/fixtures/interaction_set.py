import numpy as np
import pytest

from compgenlnc.consts import IS_POSITIVE_FIELD, PAIR_FIELDS
from compgenlnc.fileman import load_interactions
from compgenlnc.typing import dtype_interaction_tuple

from constants import INTERACTIONS_CSV


@pytest.fixture(
    params=[
        np.array(
            [
                ("NONHSAT000002.2", "hsa-let-7a-1", True),
                ("NONHSAT000002.2", "hsa-let-7a-3", True),
                ("NONHSAT000003.2", "hsa-let-7e-5p", True),
                ("NONHSAT000004.2", "hsa-let-7e", True),
                ("NONHSAT000005.2", "hsa-let-7a-1", True),
                ("NONHSAT000007.2", "hsa-let-7g", True),
            ],
            dtype_interaction_tuple,
        ),
        np.array(
            [
                ("NONHSAT000002.2", "hsa-let-7e-3p", False),
                ("NONHSAT000003.2", "hsa-let-7a-2-3p", False),
                ("NONHSAT000005.2", "hsa-let-7b-5p", False),
                ("NONHSAT000006.2", "hsa-let-7f-5p", False),
                ("NONHSAT000007.2", "hsa-let-7d-3p", False),
            ],
            dtype_interaction_tuple,
        ),
        np.array(
            [
                ("NONHSAT000005.2", "hsa-let-7a-1", True),
                ("NONHSAT000007.2", "hsa-let-7g", True),
                ("NONHSAT000008.2", "hsa-let-7c", True),
                ("NONHSAT000009.2", "hsa-let-7a-5p", True),
                ("NONHSAT000010.2", "hsa-let-7e-3p", True),
                ("NONHSAT000011.2", "hsa-let-7b-5p", True),
                ("NONHSAT000002.2", "hsa-let-7e-3p", False),
                ("NONHSAT000003.2", "hsa-let-7a-2-3p", False),
                ("NONHSAT000005.2", "hsa-let-7b-5p", False),
                ("NONHSAT000006.2", "hsa-let-7f-5p", False),
                ("NONHSAT000007.2", "hsa-let-7d-3p", False),
            ],
            dtype_interaction_tuple,
        ),
    ]
)
def interaction_set(request):
    return request.param


@pytest.fixture
def interaction_list():
    return load_interactions(INTERACTIONS_CSV)


@pytest.fixture
def interaction_positive_pairs(interaction_set):
    mask = interaction_set[IS_POSITIVE_FIELD]
    return interaction_set[mask][PAIR_FIELDS]


@pytest.fixture
def interaction_negative_pairs(interaction_set):
    mask = ~interaction_set[IS_POSITIVE_FIELD]
    return interaction_set[mask][PAIR_FIELDS]
