import pytest

from compgenlnc.extractors.molecules import (
    iter_seq_fasta
)
from constants import SEQ_FASTA_EXAMPLE


@pytest.fixture
def seq_fasta_iterator():
    return iter_seq_fasta(SEQ_FASTA_EXAMPLE)


@pytest.fixture
def seq_fasta_list(scope='module'):
    return list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))