import pytest

from compgenlnc.utils import iter_seq_ss, load_fasta

from constants import (
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PRE_MIRNA_FASTA,
    SEQ_SS_FOLDER,
)
from params import id_list


@pytest.fixture(params=list(iter_seq_ss(SEQ_SS_FOLDER)))
def seq_ss_record(request):
    return request.param


@pytest.fixture(scope="session")
def seq_ss_list():
    return list(iter_seq_ss(SEQ_SS_FOLDER))


@pytest.fixture(
    params=[
        list(iter_seq_ss(SEQ_SS_FOLDER)),
        list(iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [0, 1, 2]])),
        list(
            iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [2, 4, 6, 8]])
        ),
        list(iter_seq_ss(SEQ_SS_FOLDER, keep=[id_list[i] for i in [7, 8]])),
    ]
)
def seq_ss_filtered_list(request):
    return request.param


@pytest.fixture(scope="session")
def lnc_records():
    return list(load_fasta(LNCRNA_FASTA))


@pytest.fixture(scope="session")
def pre_mir_records():
    return list(load_fasta(PRE_MIRNA_FASTA))


@pytest.fixture(scope="session")
def mat_mir_records():
    return list(load_fasta(MAT_MIRNA_FASTA))


@pytest.fixture(scope="session")
def mir_records(pre_mir_records, mat_mir_records):
    return pre_mir_records + mat_mir_records
