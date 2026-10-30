import pytest

import compgenlnc as cgl


@pytest.fixture(scope="session")
def ss_record_identified_list():
    return list(cgl.load_fasta(cgl.paths.MIRNA_SS_FASTA, mode="ss"))


@pytest.fixture(scope="session")
def ss_fasta_list():
    return list(cgl.load_fasta(cgl.paths.MIRNA_SS_FASTA, mode="ss"))
