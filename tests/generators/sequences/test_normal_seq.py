import math

import numpy as np
import pytest

from compgenlnc.generators.sequences.normal_sequence import (
    gen_normalized_sequence_dict,
    gen_normalized_sequence_dict_from_folder,
    normalize_sequence_list,
    recode_sequence,
    recode_sequence_list,
)
from compgenlnc.structs import SeqRecord, SeqSSRecord
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import iter_seq_ss

from constants import SEQ_SS_FOLDER


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            "AACCGGTT",
            np.array([1, 1, 2, 2, 3, 3, 4, 4], np.uint8),
        ),
        (
            SeqRecord("UACGUGCAU"),
            np.array([4, 1, 2, 3, 4, 3, 2, 1, 4], np.uint8),
        ),
        (
            SeqSSRecord("", "CGCGCGCG", ""),
            np.array([2, 3] * 4, np.uint8),
        ),
    ],
)
def test_recode_sequence(seq, expected):
    recoded = recode_sequence(seq)
    assert np.array_equal(recoded, expected)


def test_recode_sequence_list(seq_fasta_list):
    recoded_list = recode_sequence_list(seq_fasta_list)
    for i, (record, expected) in enumerate(zip(seq_fasta_list, recoded_list)):
        assert np.array_equal(recode_sequence(record.seq), expected)


@pytest.mark.parametrize("normal_size", [None, 100, 50, 73])
def test_normalize_sequence_list(normal_size, seq_fasta_list):
    recoded_list = list(recode_sequence_list(seq_fasta_list))
    normal_matrix = normalize_sequence_list(recoded_list, normal_size)
    if normal_size is None:
        normal_size = math.ceil(np.mean([len(seq) for seq in recoded_list]))
    assert normal_matrix.shape == (len(seq_fasta_list), normal_size)

    average = abs(normal_matrix.mean(0))
    deviation = abs(normal_matrix.std(0))
    assert np.all(average < 1e-6)
    assert np.all((abs(deviation - 1) < 1e-6) | (deviation < 1e-10))


def test_gen_normalized_sequence_dict(seq_ss_filtered_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_normalized_sequence_dict(seq_ss_filtered_list, dict_file)
    assert dict_file.exists()

    recoded_list = list(recode_sequence_list(seq_ss_filtered_list))
    normal_matrix = normalize_sequence_list(recoded_list)
    normal_dict = load_dict(dict_file)

    for i, key in enumerate(normal_dict.keys()):
        assert np.array_equal(normal_matrix[i], normal_dict[key])


def test_gen_normalized_sequence_dict_from_folder(mirna_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_normalized_sequence_dict_from_folder(
        SEQ_SS_FOLDER, dict_file, keep=mirna_id_filter
    )
    assert dict_file.exists()

    recoded_list = recode_sequence_list(
        iter_seq_ss(SEQ_SS_FOLDER, keep=mirna_id_filter)
    )
    normal_matrix = normalize_sequence_list(recoded_list)
    normal_dict = load_dict(dict_file)

    for i, key in enumerate(normal_dict.keys()):
        assert np.array_equal(normal_matrix[i], normal_dict[key])
