import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def ss_record_identified_list():
    return cgl.load_fasta(cgl.paths.MIRNA_SS_FASTA, mode="ss")


@pytest.fixture(scope="module")
def ss_fasta_list():
    return cgl.load_fasta(cgl.paths.MIRNA_SS_FASTA, mode="ss")
