import pytest

import compgenlnc as cgl


@pytest.fixture(scope="session")
def index_list():
    return cgl.load_indexes()
