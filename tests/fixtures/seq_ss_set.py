import pytest

import compgenlnc as cgl


@pytest.fixture(scope="module")
def seq_ss_list():
    return cgl.load_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)


@pytest.fixture(scope="module")
def lnc_records():
    return cgl.load_seq_ss(cgl.paths.LNCRNA_SEQ_SS_DIR)


@pytest.fixture(scope="module")
def pre_mir_records():
    return cgl.load_seq_ss(cgl.paths.MIRNA_SEQ_SS_DIR)


@pytest.fixture(scope="module")
def mat_mir_records():
    return cgl.load_fasta(cgl.paths.MAT_MIRNA_FASTA)


@pytest.fixture(scope="module")
def mir_records(mat_mir_records):
    return mat_mir_records
