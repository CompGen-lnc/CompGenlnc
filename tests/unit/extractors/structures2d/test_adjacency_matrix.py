import math

import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "ss, expected",
    [
        (
            "(())",
            np.array(
                [
                    [0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0],
                ],
                np.bool,
            ),
        ),
        (
            cgl.structs.StructureRecord("(..)."),
            np.array(
                [
                    [0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                ],
                np.bool,
            ),
        ),
        (
            cgl.structs.SeqSSRecord("", "", "(.())."),
            np.array(
                [
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ],
                np.bool,
            ),
        ),
    ],
)
def test_get_adjacency_matrix(ss, expected):
    matrix = cgl.get_adjacency_matrix(ss)
    assert np.array_equal(matrix, expected)


def test_get_adjacency_matrix_list(ss_fasta_list):
    matrix_it = cgl.get_adjacency_matrix_list(ss_fasta_list)
    for i, (record, expected) in enumerate(zip(ss_fasta_list, matrix_it)):
        assert np.array_equal(cgl.get_adjacency_matrix(record.ss), expected)


@pytest.mark.parametrize("normal_size", [None, 100, 50, 73])
def test_normalize_matrix_list(normal_size, ss_fasta_list):
    matrix_it = cgl.get_adjacency_matrix_list(ss_fasta_list)
    normal_it = cgl.normalize_matrix_list(matrix_it, normal_size)
    if normal_size is None:
        normal_size = math.ceil(np.mean([len(ss.ss) for ss in ss_fasta_list]))

    for mat in normal_it:
        assert mat.shape == (normal_size, normal_size)


def test_gen_normalized_matrix_dict(seq_ss_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_normalized_matrix_dict(seq_ss_list, dict_file)
    assert dict_file.exists()

    matrix_it = cgl.get_adjacency_matrix_list(seq_ss_list)
    normal_it = cgl.normalize_matrix_list(matrix_it)
    normal_dict = cgl.load_dict(dict_file)

    for mat, expected in zip(normal_dict.values(), normal_it):
        assert np.array_equal(mat, expected)


def test_gen_normalized_matrix_dict_from_folder(example_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_normalized_matrix_dict_from_folder(
        cgl.paths.MIRNA_SEQ_SS_FOLDER, dict_file, keep=example_id_filter
    )
    assert dict_file.exists()

    matrix_it = cgl.get_adjacency_matrix_list(
        cgl.iter_seq_ss(cgl.paths.MIRNA_SEQ_SS_FOLDER, keep=example_id_filter)
    )
    normal_it = cgl.normalize_matrix_list(matrix_it)
    normal_dict = cgl.load_dict(dict_file)

    for mat, expected in zip(normal_dict.values(), normal_it):
        assert np.array_equal(mat, expected)
