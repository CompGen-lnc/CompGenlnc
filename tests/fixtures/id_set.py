import pytest

from compgenlnc.utils import load_fasta

from constants import (
    SEQ_FASTA_EXAMPLE,
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PRE_MIRNA_FASTA,
)


@pytest.fixture(params=[rec.id for rec in load_fasta(SEQ_FASTA_EXAMPLE)])
def example_id(request):
    return request.param


@pytest.fixture(
    params=[
        [],
        [
            rec.id
            for i, rec in enumerate(load_fasta(SEQ_FASTA_EXAMPLE))
            if i in [0, 1, 2]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(SEQ_FASTA_EXAMPLE))
            if i in [2, 4, 6, 8]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(SEQ_FASTA_EXAMPLE))
            if i in [7, 8]
        ],
    ]
)
def example_id_filter(request):
    return request.param


@pytest.fixture(scope="session")
def lncrna_id_list():
    records = load_fasta(LNCRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(
    params=[
        [],
        [
            rec.id
            for i, rec in enumerate(load_fasta(LNCRNA_FASTA))
            if i in [0, 1, 2]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(LNCRNA_FASTA))
            if i in [2, 4, 6, 8]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(LNCRNA_FASTA))
            if i in [7, 8]
        ],
    ],
)
def lncrna_id_filter(request):
    return request.param


@pytest.fixture(scope="session")
def pre_mirna_id_list():
    records = load_fasta(PRE_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(
    params=[
        [],
        [
            rec.id
            for i, rec in enumerate(load_fasta(PRE_MIRNA_FASTA))
            if i in [0, 1, 2]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(PRE_MIRNA_FASTA))
            if i in [2, 4, 6, 8]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(PRE_MIRNA_FASTA))
            if i in [7, 8]
        ],
    ],
)
def pre_mirna_id_filter(request):
    return request.param


@pytest.fixture(scope="session")
def mat_mirna_id_list():
    records = load_fasta(MAT_MIRNA_FASTA)
    return [rec.id for rec in records]


@pytest.fixture(
    params=[
        [],
        [
            rec.id
            for i, rec in enumerate(load_fasta(MAT_MIRNA_FASTA))
            if i in [0, 1, 2, 10, 11]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(MAT_MIRNA_FASTA))
            if i in [2, 4, 6, 8, 10, 12, 14]
        ],
        [
            rec.id
            for i, rec in enumerate(load_fasta(MAT_MIRNA_FASTA))
            if i in [7, 8, 13, 14]
        ],
    ],
)
def mat_mirna_id_filter(request):
    return request.param


@pytest.fixture(scope="session")
def mirna_id_list(pre_mirna_id_list, mat_mirna_id_list):
    return pre_mirna_id_list + mat_mirna_id_list


@pytest.fixture
def mirna_id_filter(pre_mirna_id_filter, mat_mirna_id_filter):
    return pre_mirna_id_filter + mat_mirna_id_filter
