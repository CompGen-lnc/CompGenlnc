import pytest

from compgenlnc.utils.fasta_manager import (
    load_fasta,
    iter_seq_fasta,
    iter_seq_ss
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.fixture
def seq_fasta_loader():
    return list(load_fasta(SEQ_FASTA_EXAMPLE, mode='seq'))


@pytest.fixture
def seq_ss_list(scope='session'):
    return list(iter_seq_ss(SEQ_SS_FOLDER))