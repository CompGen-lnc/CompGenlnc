import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss,
    read_sequence
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.mark.parametrize(
    'expected, index',
    [
        (SeqSSRecord('hsa1', '', ''), 0),
        (SeqSSRecord('hsa2', 'ACGUGCA', ''), 1),
        (SeqSSRecord('mol1', 'ACGUACGUACGU', ''), 2),
        (SeqSSRecord('mol2', 'AAAA', ''), 3),
        (SeqSSRecord('mol3', 'GCGCGCGC', ''), 4),
    ]
)
def test_iter_seq_fasta(expected, index, seq_fasta_list):
    record = seq_fasta_list[index]
    assert record.id == expected.id
    assert record.seq == expected.seq


@pytest.mark.parametrize(
    'filter, expected',
    [
        (['mol'], 3),
        (['hsa'], 2),
        (['mol', 'hsa'], 5),
        (['non'], 0),
    ]
)
def test_iter_seq_fasta_filtered(filter, expected):
    records_list = [
        record
        for record in iter_seq_fasta(SEQ_FASTA_EXAMPLE, filter)
    ]
    assert len(records_list) == expected


@pytest.mark.parametrize(
    'expected, index',
    [
        (SeqSSRecord('hsa1', '', ''), 0),
        (SeqSSRecord('hsa2', 'ACGTGCA', ''), 1),
        (SeqSSRecord('mol1', 'ACGTACGTACGT', ''), 2),
        (SeqSSRecord('mol2', 'AAAA', ''), 3),
        (SeqSSRecord('mol3', 'GCGCGCGC', ''), 4),
    ]
)
def test_iter_seq_ss(expected, index, seq_ss_list):
    record = seq_ss_list[index]
    assert record.id == expected.id
    assert record.seq == expected.seq


@pytest.mark.parametrize(
    'keep, expected',
    [
        (['mol1', 'hsa2', 'mol2'], 3),
        (['hsa1', 'hsa2'], 2),
        ([], 5),
        (['non1'], 0),
    ]
)
def test_iter_seq_ss_filtered(keep, expected):
    records_list = [
        record
        for record in iter_seq_ss(SEQ_SS_FOLDER, keep=keep)
    ]
    assert len(records_list) == expected


@pytest.mark.parametrize(
    'id_, expected',
    [
        ('hsa1', ''),
        ('hsa2', 'ACGTGCA'),
        ('mol1', 'ACGTACGTACGT'),
        ('mol2', 'AAAA'),
        ('mol3', 'GCGCGCGC'),
    ]
)
def test_read_sequence(id_, expected):
    filename = SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{id_}.dat'
    assert filename.exists()

    seq = read_sequence(filename)
    assert seq == expected
