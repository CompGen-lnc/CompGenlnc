import pytest

from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.fixture
def seq_fasta_iterator():
    return iter_seq_fasta(SEQ_FASTA_EXAMPLE)


@pytest.fixture
def seq_fasta_list(scope='session'):
    return list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))


@pytest.fixture
def seq_ss_iterator():
    return iter_seq_ss(SEQ_SS_FOLDER)


@pytest.fixture
def seq_ss_list(scope='session'):
    return list(iter_seq_ss(SEQ_SS_FOLDER))