import os

import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.extractors.molecules.fasta_constructor import iter_seq_fasta
from compgenlnc.extractors.molecules.seq_ss_constructor import (
    fasta_to_seq_ss,
    join_seq_ss,
)
from constants import (
    SEQ_FASTA_EXAMPLE
)


def test_fasta_to_seq_ss(tmp_path):
    folder = tmp_path / 'seq+ss_join'
    total = fasta_to_seq_ss(SEQ_FASTA_EXAMPLE, folder)
    assert total == 5


@pytest.mark.parametrize(
    'record',
    list(iter_seq_fasta(SEQ_FASTA_EXAMPLE))
)
def test_join_seq_ss(record, tmp_path): 
    seq = str(record.seq)
    ss = ''
    id_ = record.id
    folder = tmp_path / 'seq+ss_join'
    folder.mkdir(parents=True, exist_ok=True)
    join_seq_ss(seq, ss, id_, folder)

    path = folder / f'{SEQ_SS_PREFIX}{id_}.txt'
    assert path.exists()

    text = path.read_text()
    assert text == f'{seq}\n{ss}'
