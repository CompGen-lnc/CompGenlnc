from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.utils.fasta_manager import (
    iter_seq_fasta,
    iter_seq_ss,
    read_sequence
)
from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.mark.parametrize(
    'expected, index',
    [
        (SeqRecord(Seq(''), 'hsa1', description=''), 0),
        (SeqRecord(Seq('ACGUGCA'), 'hsa2', description=''), 1),
        (SeqRecord(Seq('ACGUACGUACGU'), 'mol1', description=''), 2),
        (SeqRecord(Seq('AAAA'), 'mol2', description=''), 3),
        (SeqRecord(Seq('GCGCGCGC'), 'mol3', description=''), 4),
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
        (SeqRecord(Seq(''), 'hsa1', description=''), 0),
        (SeqRecord(Seq('ACGTGCA'), 'hsa2', description=''), 1),
        (SeqRecord(Seq('ACGTACGTACGT'), 'mol1', description=''), 2),
        (SeqRecord(Seq('AAAA'), 'mol2', description=''), 3),
        (SeqRecord(Seq('GCGCGCGC'), 'mol3', description=''), 4),
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
