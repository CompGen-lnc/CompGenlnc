import numpy as np
import pytest
from constants import SEQ_SS_FOLDER

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.generators.sequences.kmer import (
    gen_kmer_dict,
    gen_kmer_dict_from_folder,
    get_kmer,
    get_kmer_by_name,
)
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import read_sequence


@pytest.mark.parametrize(
    "seq, k, expected",
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
            np.array([0.25] * 4, np.float32),
        ),
    ],
)
def test_get_kmer(seq, k, expected):
    kmers = get_kmer(seq, k)
    assert np.array_equal(kmers, expected)


def test_get_kmer_by_name(mirna_id):
    kmers = get_kmer_by_name(mirna_id, SEQ_SS_FOLDER, 3)
    seq = read_sequence(SEQ_SS_FOLDER / f"{SEQ_SS_PREFIX}{mirna_id}.dat")
    expected = get_kmer(seq, 3)
    assert np.array_equal(kmers, expected)


def test_gen_kmer_dict(seq_ss_filtered_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_kmer_dict(seq_ss_filtered_list, dict_file, 3)
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all(
        [
            np.array_equal(kmer_dict[record.id], get_kmer(record.seq, 3))
            for record in seq_ss_filtered_list
        ]
    )


def test_gen_kmer_dict_from_folder(mirna_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_kmer_dict_from_folder(
        SEQ_SS_FOLDER, dict_file, 3, keep=mirna_id_filter
    )
    assert dict_file.exists()

    kmer_dict = load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                kmer_dict[id_], get_kmer_by_name(id_, SEQ_SS_FOLDER, 3)
            )
            for id_ in mirna_id_filter
        ]
    )
