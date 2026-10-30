import pytest

import compgenlnc as cgl


@pytest.fixture(scope="session")
def seq_ss_list():
    return list(cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_FOLDER))


@pytest.fixture(scope="session")
def lnc_records():
    return list(cgl.load_fasta(cgl.paths.LNCRNA_FASTA))


@pytest.fixture(scope="session")
def pre_mir_records():
    return list(cgl.load_fasta(cgl.paths.PRE_MIRNA_FASTA))


@pytest.fixture(scope="session")
def mat_mir_records():
    return list(cgl.load_fasta(cgl.paths.MAT_MIRNA_FASTA))


@pytest.fixture(scope="session")
def mir_records(mat_mir_records):
    return mat_mir_records
