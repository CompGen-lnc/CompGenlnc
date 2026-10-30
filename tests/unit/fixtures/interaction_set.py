import random

import pytest

import compgenlnc as cgl


@pytest.fixture(scope="session")
def interaction_list():
    return cgl.load_interactions(cgl.paths.PAIRS_FILE)


@pytest.fixture(
    params=[
        sorted(random.sample(range(300), 100)),
        sorted(random.sample(range(300, 600), 100)),
        sorted(random.sample(range(600), 250)),
    ]
)
def interaction_set(request, interaction_list):
    return interaction_list[request.param]


@pytest.fixture
def interaction_positive_pairs(interaction_set):
    mask = interaction_set[cgl.IS_POSITIVE_FIELD]
    return interaction_set[mask][cgl.PAIR_FIELDS]


@pytest.fixture
def interaction_negative_pairs(interaction_set):
    mask = ~interaction_set[cgl.IS_POSITIVE_FIELD]
    return interaction_set[mask][cgl.PAIR_FIELDS]
