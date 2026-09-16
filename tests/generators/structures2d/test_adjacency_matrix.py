import math

import numpy as np
import pytest

from compgenlnc.generators.structures2d import (
    normalize_matrix_list,
    gen_normalized_matrix_dict,
    gen_normalized_matrix_dict_from_fasta,
    gen_normalized_matrix_dict_from_folder,
    get_adjacency_matrix,
    get_adjacency_matrix_list,
)
from compgenlnc.structs import SeqSSRecord, StructureRecord
from compgenlnc.utils import iter_seq_ss, load_dict

from constants import SEQ_SS_FOLDER


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
            StructureRecord("(..)."),
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
            SeqSSRecord("", "", "(.())."),
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
    matrix = get_adjacency_matrix(ss)
    assert np.array_equal(matrix, expected)


def test_get_adjacency_matrix_list(ss_fasta_list):
    matrix_it = get_adjacency_matrix_list(ss_fasta_list)
    for i, (record, expected) in enumerate(zip(ss_fasta_list, matrix_it)):
        assert np.array_equal(get_adjacency_matrix(record.ss), expected)


@pytest.mark.parametrize("normal_size", [None, 100, 50, 73])
def test_normalize_matrix_list(normal_size, ss_fasta_list):
    matrix_it = get_adjacency_matrix_list(ss_fasta_list)
    normal_it = normalize_matrix_list(matrix_it, normal_size)
    if normal_size is None:
        normal_size = math.ceil(np.mean([len(ss.ss) for ss in ss_fasta_list]))

    for mat in normal_it:
        assert mat.shape == (normal_size, normal_size)


def test_gen_normalized_matrix_dict(seq_ss_filtered_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_normalized_matrix_dict(seq_ss_filtered_list, dict_file)
    assert dict_file.exists()

    matrix_it = get_adjacency_matrix_list(seq_ss_filtered_list)
    normal_it = normalize_matrix_list(matrix_it)
    normal_dict = load_dict(dict_file)

    for mat, expected in zip(normal_dict.values(), normal_it):
        assert np.array_equal(mat, expected)


def test_gen_normalized_matrix_dict_from_folder(mirna_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    gen_normalized_matrix_dict_from_folder(
        SEQ_SS_FOLDER, dict_file, keep=mirna_id_filter
    )
    assert dict_file.exists()

    matrix_it = get_adjacency_matrix_list(
        iter_seq_ss(SEQ_SS_FOLDER, keep=mirna_id_filter)
    )
    normal_it = normalize_matrix_list(matrix_it)
    normal_dict = load_dict(dict_file)

    for mat, expected in zip(normal_dict.values(), normal_it):
        assert np.array_equal(mat, expected)
