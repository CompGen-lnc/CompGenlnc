import pytest

from compgenlnc.utils import load_fasta

from constants import LNCRNA_FASTA, MAT_MIRNA_FASTA, PRE_MIRNA_FASTA
from params import id_list


@pytest.fixture(params=id_list)
def example_id(request):
    return request.param


@pytest.fixture(
    params=[
        [],
        [id_list[i] for i in [0, 1, 2]],
        [id_list[i] for i in [2, 4, 6, 8]],
        [id_list[i] for i in [7, 8]],
    ]
)
def example_id_filter(request):
    return request.param


@pytest.fixture(scope="session")
def lncrna_id_list():
    records = load_fasta(LNCRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(scope="session")
def pre_mirna_id_list():
    records = load_fasta(PRE_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(scope="session")
def mat_mirna_id_list():
    records = load_fasta(MAT_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(scope="session")
def mirna_id_list(pre_mirna_id_list, mat_mirna_id_list):
    return pre_mirna_id_list + mat_mirna_id_list
