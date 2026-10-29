import pytest
from constants import EMBEDDING_FOLDER

from compgenlnc import load_indexes


@pytest.fixture(scope="session")
def index_list():
    return load_indexes(index_file=EMBEDDING_FOLDER / "index.npy")
