from pathlib import Path

import pytest

import compgenlnc as cgl

pytestmark = pytest.mark.unit


@pytest.fixture(scope="package", autouse=True)
def set_base_path():
    original_path = cgl.config.BASE_PATH
    unit_path = Path(__file__).resolve().parent
    cgl.set_base_path(unit_path)
    cgl.dict_indexer._cgl_indexes = None
    cgl.dict_indexer._lnc_indexed_records = None
    cgl.dict_indexer._mat_mir_indexed_records = None
    cgl.dict_indexer._pre_mir_indexed_records = None

    yield

    cgl.set_base_path(original_path)
