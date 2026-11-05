import math

import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            "AACCGGTT",
            np.array([1, 1, 2, 2, 3, 3, 4, 4], np.uint8),
        ),
        (
            cgl.structs.SeqRecord("UACGUGCAU"),
            np.array([4, 1, 2, 3, 4, 3, 2, 1, 4], np.uint8),
        ),
        (
            cgl.structs.SeqSSRecord("", "CGCGCGCG", ""),
            np.array([2, 3] * 4, np.uint8),
        ),
    ],
)
def test_recode_sequence(seq, expected):
    recoded = cgl.recode_sequence(seq)
    assert np.array_equal(recoded, expected)


def test_recode_sequence_list(seq_fasta_list):
    recoded_list = cgl.recode_sequence_list(seq_fasta_list)
    for i, (record, expected) in enumerate(zip(seq_fasta_list, recoded_list)):
        assert np.array_equal(cgl.recode_sequence(record.seq), expected)


@pytest.mark.parametrize("normal_size", [None, 100, 50, 73])
def test_normalize_sequence_list(normal_size, seq_fasta_list):
    recoded_list = list(cgl.recode_sequence_list(seq_fasta_list))
    normal_matrix = cgl.normalize_sequence_list(recoded_list, normal_size)
    if normal_size is None:
        normal_size = math.ceil(np.mean([len(seq) for seq in recoded_list]))
    assert normal_matrix.shape == (len(seq_fasta_list), normal_size)

    average = abs(normal_matrix.mean(0))
    deviation = abs(normal_matrix.std(0))
    assert np.all(average < 1e-6)
    assert np.all((abs(deviation - 1) < 1e-5) | (deviation < 1e-10))


def test_gen_normalized_sequence_dict(seq_ss_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_normalized_sequence_dict(dict_file, seq_ss_list)
    assert dict_file.exists()

    recoded_list = list(cgl.recode_sequence_list(seq_ss_list))
    normal_matrix = cgl.normalize_sequence_list(recoded_list)
    normal_dict = cgl.load_dict(dict_file)

    for i, key in enumerate(normal_dict.keys()):
        assert np.array_equal(normal_matrix[i], normal_dict[key])


def test_gen_normalized_sequence_dict_from_folder(example_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_normalized_sequence_dict_from_folder(
        dict_file, cgl.paths.MIRNA_SEQ_SS_FOLDER, keep=example_id_filter
    )
    assert dict_file.exists()

    recoded_list = cgl.recode_sequence_list(
        cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_FOLDER, keep=example_id_filter)
    )
    normal_matrix = cgl.normalize_sequence_list(recoded_list)
    normal_dict = cgl.load_dict(dict_file)

    for i, key in enumerate(normal_dict.keys()):
        assert np.array_equal(normal_matrix[i], normal_dict[key])
