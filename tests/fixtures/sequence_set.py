import pytest
from constants import SEQ_FASTA_EXAMPLE

from compgenlnc.utils.fasta_manager import iter_seq_fasta, load_fasta


@pytest.fixture(scope="session")
def seq_fasta_loader():
    return list(load_fasta(SEQ_FASTA_EXAMPLE, mode="seq"))


@pytest.fixture(scope="session")
def seq_fasta_list():
    return list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))
