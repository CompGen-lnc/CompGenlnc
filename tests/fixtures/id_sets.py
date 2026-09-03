import pytest

from params import id_list


@pytest.fixture(
    params=id_list
)
def mirna_id(request):
    return request.param


@pytest.fixture(
    params=[
        [],
        [id_list[i] for i in [0, 1, 2]],
        [id_list[i] for i in [2, 4, 6, 8]],
        [id_list[i] for i in [7, 8]],
    ]
)
def mirna_id_filter(request):
    return request.param