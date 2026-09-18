import pytest
from constants import SS_FASTA_EXAMPLE

from compgenlnc.utils.fasta_manager import load_fasta


@pytest.fixture(
    params=list(load_fasta(SS_FASTA_EXAMPLE, mode="ss")),
)
def ss_record_identified(request):
    return request.param


@pytest.fixture(scope="session")
def ss_fasta_list():
    return list(load_fasta(SS_FASTA_EXAMPLE, mode="ss"))
