import pytest

from compgenlnc.utils.fasta_manager import load_fasta, iter_seq_fasta

from constants import SEQ_FASTA_EXAMPLE


@pytest.fixture(scope='session')
def seq_fasta_loader():
    return list(load_fasta(SEQ_FASTA_EXAMPLE, mode='seq'))


@pytest.fixture(scope='session')
def seq_fasta_list():
    return list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))
