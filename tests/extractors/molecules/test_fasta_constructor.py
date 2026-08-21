from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from compgenlnc.extractors.molecules import filter_fasta, iter_seq_fasta
from constants import SEQ_FASTA_EXAMPLE


@pytest.mark.parametrize(
    'expected, index',
    [
        (SeqRecord(Seq('ACGUACGUACGU'), 'mol1', description=''), 0),
        (SeqRecord(Seq('AAAA'), 'mol2', description=''), 1),
        (SeqRecord(Seq('GCGCGCGC'), 'mol3', description=''), 2),
        (SeqRecord(Seq(''), 'hsa1', description=''), 3),
        (SeqRecord(Seq('ACGUGCA'), 'hsa2', description=''), 4),
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
    'filter, expected',
    [
        (['mol'], ['mol1', 'mol2', 'mol3']),
        (['hsa'], ['hsa1', 'hsa2']),
        (['mol', 'hsa'], ['hsa1', 'hsa2', 'mol1', 'mol2', 'mol3']),
        (['non'], []),
    ]
)
def test_filter_seq_fasta(filter, expected, tmp_path):
    new_fasta = tmp_path / 'new_fasta.fa'
    filter_fasta(SEQ_FASTA_EXAMPLE, filter, new_fasta)
    assert new_fasta.exists()

    fasta = new_fasta.read_text()
    ids = [record.id for record in iter_seq_fasta(new_fasta)]
    assert ids == expected