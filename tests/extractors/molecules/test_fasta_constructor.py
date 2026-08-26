import pytest

from compgenlnc.extractors.molecules.fasta_constructor import (
    filter_fasta,
    seq_ss_to_fasta
)
from compgenlnc.utils.fasta_manager import iter_seq_fasta

from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER


@pytest.mark.parametrize(
    'filter, expected',
    [
        (['hsa'], ['hsa1', 'hsa2']),
        (['mol'], ['mol1', 'mol2', 'mol3']),
        (['mol', 'hsa'], ['hsa1', 'hsa2', 'mol1', 'mol2', 'mol3']),
        (['non'], []),
    ]
)
def test_filter_seq_fasta(filter, expected, tmp_path):
    new_fasta = tmp_path / 'new_fasta.fa'
    filter_fasta(SEQ_FASTA_EXAMPLE, new_fasta, filter)
    assert new_fasta.exists()

    ids = [record.id for record in iter_seq_fasta(new_fasta)]
    assert ids == expected


def test_seq_ss_to_fasta(tmp_path):
    new_seq_fasta_file = tmp_path / 'new_seq_fasta.fa'
    new_ss_fasta_file = tmp_path / 'new_ss_fasta.fa'
    seq_ss_to_fasta(SEQ_SS_FOLDER, new_seq_fasta_file, new_ss_fasta_file)
    assert new_seq_fasta_file.exists()

    new_seq_fasta = new_seq_fasta_file.read_text()
    seq_fasta = SEQ_FASTA_EXAMPLE.read_text()
    assert new_seq_fasta.upper() == seq_fasta.upper().replace('U', 'T')
