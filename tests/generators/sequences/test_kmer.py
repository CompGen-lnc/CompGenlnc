import numpy as np
import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.generators.sequences.kmer import (
    gen_kmer_dict,
    gen_kmer_dict_from_folder,
    get_kmer,
    get_kmer_by_name,
)
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import iter_seq_ss, read_sequence
from constants import SEQ_SS_FOLDER

@pytest.mark.parametrize(
    'seq, k, expected',
    [
        ('ACGTGAAC', 2, np.array([
            0.375, 0.125, 0.25, 0.25, 1 / 7, 0, 2 / 7, 0, 0,
            0, 0, 1 / 7, 0, 0, 0, 1 / 7, 1 / 7, 1 / 7, 0, 0,
        ], np.float32)),
        ('AAA', 3, np.array(
            [1] + [0] * 3 + [1] + [0] * 15 + [1] + [0] * 63,
            np.float32
        )),
        (
            SeqSSRecord('', 'ATCGATCGATCG', ''),
            1,
            np.array([0.25]*4, np.float32)
        ),
    ]
)
def test_get_kmer(seq, k, expected):
    kmers = get_kmer(seq, k)
    assert np.array_equal(kmers, expected)


@pytest.mark.parametrize(
    'id_, k', [('hsa1', 3), ('hsa2', 1), ('mol1', 2), ('mol2', 3), ('mol3', 1)]
)
def test_get_kmer_by_name(id_, k):
    kmers = get_kmer_by_name(id_, SEQ_SS_FOLDER, k)
    seq = read_sequence(SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{id_}.dat')
    expected = get_kmer(seq, k)
    assert np.array_equal(kmers, expected)


@pytest.mark.parametrize(
    'record_list, k',
    [
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'hsa2', 'mol1']), 3),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'mol3', 'mol2']), 2),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa2', 'mol2']), 4),
        (iter_seq_ss(SEQ_SS_FOLDER, keep=['hsa1', 'hsa2', 'mol1']), 1),
    ]
)
def test_gen_kmer_dict(record_list, k, tmp_path):
    dict_file = tmp_path / 'kmer.dict'
    gen_kmer_dict(record_list, dict_file, k)
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            kmer_dict[record.id],
            get_kmer(record.seq, k)
        )
        for record in record_list
    ])


@pytest.mark.parametrize(
    'seq_list, k',
    [
        (['hsa1', 'hsa2', 'mol1'], 3),
        (['hsa1', 'mol3', 'mol2'], 2),
        (['hsa2', 'mol2'], 4),
        (['hsa1', 'hsa2', 'mol1'], 1),
    ]
)
def test_gen_kmer_dict_from_folder(seq_list, k, tmp_path):
    dict_file = tmp_path / 'kmer.dict'
    gen_kmer_dict_from_folder(SEQ_SS_FOLDER, dict_file, k, keep=seq_list)
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            kmer_dict[id_],
            get_kmer_by_name(id_, SEQ_SS_FOLDER, k)
        )
        for id_ in seq_list
    ])

