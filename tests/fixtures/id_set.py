import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def example_id_list():
    records = cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(params=[[], range(50), range(0, 200, 2), range(300, 350)])
def example_id_filter(request, example_id_list):
    filter = request.param
    return [example_id_list[i] for i in filter]


@pytest.fixture(scope="module")
def lncrna_id_list():
    records = cgl.load_fasta(cgl.paths.LNCRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(params=[[], range(50), range(0, 200, 2), range(300, 350)])
def lncrna_id_filter(request, lncrna_id_list):
    filter = request.param
    return [lncrna_id_list[i] for i in filter]


@pytest.fixture(scope="module")
def pre_mirna_id_list():
    records = cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(params=[[], range(50), range(0, 200, 2), range(300, 350)])
def pre_mirna_id_filter(request, pre_mirna_id_list):
    filter = request.param
    return [pre_mirna_id_list[i] for i in filter]


@pytest.fixture(scope="module")
def mat_mirna_id_list():
    records = cgl.load_fasta(cgl.paths.MAT_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(params=[[], range(50), range(0, 200, 2), range(300, 350)])
def mat_mirna_id_filter(request, mat_mirna_id_list):
    filter = request.param
    return [mat_mirna_id_list[i] for i in filter]


@pytest.fixture(scope="module")
def mirna_id_list(mat_mirna_id_list):
    return mat_mirna_id_list


@pytest.fixture
def mirna_id_filter(mat_mirna_id_filter):
    return mat_mirna_id_filter


@pytest.fixture(scope="module")
def precursors_consverter():
    return cgl.load_precursors(cgl.paths.DATA_MIRNA)
