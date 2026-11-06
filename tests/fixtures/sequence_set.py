import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def seq_fasta_loader():
    return cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq")


@pytest.fixture(scope="module")
def seq_fasta_list():
    return cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA, mode="seq")
