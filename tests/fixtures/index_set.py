import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def index_list():
    return cgl.load_indexes()
