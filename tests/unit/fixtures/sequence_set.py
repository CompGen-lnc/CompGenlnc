import pytest

import compgenlnc as cgl


@pytest.fixture(scope="session")
def seq_fasta_loader():
    return list(cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq"))


@pytest.fixture(scope="session")
def seq_fasta_list():
    return list(cgl.iter_seq_fasta(cgl.paths.PRE_MIRNA_FASTA))
