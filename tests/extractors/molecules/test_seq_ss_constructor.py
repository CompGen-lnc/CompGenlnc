import os

import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.utils.fasta_manager import iter_seq_fasta, read_sequence
from compgenlnc.extractors.molecules.seq_ss_constructor import (
    fasta_to_seq_ss,
    join_seq_ss,
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER, SS_FASTA_EXAMPLE


def test_fasta_to_seq_ss(tmp_path):
    folder = tmp_path / 'seq+ss_join'
    total = fasta_to_seq_ss(SEQ_FASTA_EXAMPLE, SS_FASTA_EXAMPLE, folder)
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

    path = folder / f'{SEQ_SS_PREFIX}{id_}.dat'
    assert path.exists()

    new_seq = read_sequence(path)
    expected = read_sequence(SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{id_}.dat')
    assert new_seq == expected
