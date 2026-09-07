import pytest

from compgenlnc.utils.fasta_manager import iter_seq_ss

from constants import SEQ_SS_FOLDER
from params import id_list


@pytest.fixture(
    scope='session',
    params=list(iter_seq_ss(SEQ_SS_FOLDER)),
)
def seq_ss_record(request):
    return request.param


@pytest.fixture(scope='session')
def seq_ss_list():
    return list(iter_seq_ss(SEQ_SS_FOLDER))


@pytest.fixture(
    params=[
        list(iter_seq_ss(SEQ_SS_FOLDER)),
        list(iter_seq_ss(
            SEQ_SS_FOLDER,
            keep=[id_list[i] for i in [0, 1, 2]]
        )),
        list(iter_seq_ss(
            SEQ_SS_FOLDER,
            keep=[id_list[i] for i in [2, 4, 6, 8]]
        )),
        list(iter_seq_ss(
            SEQ_SS_FOLDER,
            keep=[id_list[i] for i in [7, 8]]
        )),
    ]
)
def seq_ss_filtered_list(request):
    return request.param