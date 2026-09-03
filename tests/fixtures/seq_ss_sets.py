import pytest

from compgenlnc.utils.fasta_manager import iter_seq_ss

from constants import SEQ_SS_FOLDER
from params import id_list


@pytest.fixture
def seq_ss_list(scope='session'):
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