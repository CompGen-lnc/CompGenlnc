import pytest

from compgenlnc import load_indexes

from constants import EMBEDDING_FOLDER


@pytest.fixture(scope="session")
def index_list():
    return load_indexes(index_file=EMBEDDING_FOLDER / "index.npy")
