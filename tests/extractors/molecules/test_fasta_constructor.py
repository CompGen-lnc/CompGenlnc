import pytest

from compgenlnc.extractors import extract_2d_structure
from compgenlnc.extractors.molecules import (
    filter_fasta,
    gen_fasta_2d,
    seq_ss_to_fasta,
)
from compgenlnc.utils import load_fasta

from constants import SEQ_FASTA_EXAMPLE, SEQ_SS_FOLDER, SS_FASTA_EXAMPLE


@pytest.mark.parametrize(
    'filter, expected',
    [
        (['7a'], ['hsa-let-7a-1', 'hsa-let-7a-2', 'hsa-let-7a-3']),
        (['7f', '7c'], ['hsa-let-7c', 'hsa-let-7f-1', 'hsa-let-7f-2']),
        (['hsa'], [record.id
                   for record in load_fasta(SEQ_FASTA_EXAMPLE, mode='seq')]),
        (['non'], []),
    ]
)
def test_filter_seq_fasta(filter, expected, tmp_path):
    new_fasta = tmp_path / 'new_fasta.fa'
    filter_fasta(SEQ_FASTA_EXAMPLE, new_fasta, filter)
    assert new_fasta.exists()

    ids = [record.id for record in load_fasta(new_fasta, mode='seq')]
    assert ids == expected


def test_gen_fasta_2d(mirna_id_filter, tmp_path):
    ss_fasta = tmp_path / 'fasta2d.fa'
    gen_fasta_2d(SEQ_FASTA_EXAMPLE, ss_fasta, keep=mirna_id_filter)
    seq_list = load_fasta(SEQ_FASTA_EXAMPLE, mode='seq', keep=mirna_id_filter)
    ss_list = load_fasta(ss_fasta, mode='ss')

    for seq, expected in zip(seq_list, ss_list):
        ss = extract_2d_structure(seq)
        assert ss == expected.ss


def test_seq_ss_to_fasta(tmp_path):
    new_seq_fasta_file = tmp_path / 'new_seq_fasta.fa'
    new_ss_fasta_file = tmp_path / 'new_ss_fasta.fa'
    seq_ss_to_fasta(SEQ_SS_FOLDER, new_seq_fasta_file, new_ss_fasta_file)
    assert new_seq_fasta_file.exists()

    new_seq_fasta = new_seq_fasta_file.read_text().replace('\n', '')
    seq_fasta = SEQ_FASTA_EXAMPLE.read_text().replace('\n', '')
    assert new_seq_fasta.upper() == seq_fasta.upper().replace('U', 'T')
