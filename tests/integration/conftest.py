from pathlib import Path

import pytest

import compgenlnc as cgl

integration_path = Path(__file__).resolve().parent


def pytest_collection_modifyitems(items):
    for item in items:
        if Path(item.fspath).is_relative_to(integration_path):
            item.add_marker("integration")


@pytest.fixture(scope="package", autouse=True)
def set_base_path():
    original_path = cgl.config.base_path
    cgl.set_base_path(integration_path)
    cgl.dict_indexer._cgl_indexes = None
    cgl.dict_indexer._lnc_indexed_records = None
    cgl.dict_indexer._mat_mir_indexed_records = None
    cgl.dict_indexer._pre_mir_indexed_records = None

    yield

    cgl.set_base_path(original_path)
