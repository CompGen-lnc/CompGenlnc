import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "seq, k, expected",
    [
        (
            "ACGTGAAC",
            2,
            np.array(
                [
                    0.375,
                    0.125,
                    0.25,
                    0.25,
                    1 / 7,
                    0,
                    2 / 7,
                    0,
                    0,
                    0,
                    0,
                    1 / 7,
                    0,
                    0,
                    0,
                    1 / 7,
                    1 / 7,
                    1 / 7,
                    0,
                    0,
                ],
                np.float32,
            ),
        ),
        (
            "AAA",
            3,
            np.array(
                [1] + [0] * 3 + [1] + [0] * 15 + [1] + [0] * 63, np.float32
            ),
        ),
        (
            cgl.structs.SeqSSRecord("", "ATCGATCGATCG", ""),
            1,
            np.array([0.25] * 4, np.float32),
        ),
    ],
)
def test_get_kmer(seq, k, expected):
    kmers = cgl.get_kmer(seq, k)
    assert np.array_equal(kmers, expected)


def test_get_kmer_by_name(example_id_list):
    for id_ in example_id_list:
        kmers = cgl.get_kmer_by_name(cgl.paths.MIRNA_SEQ_SS_FOLDER, id_, 3)
        seq = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_FOLDER, id_)
        expected = cgl.get_kmer(seq, 3)
        assert np.array_equal(kmers, expected)


def test_gen_kmer_dict(seq_ss_list, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_kmer_dict(dict_file, seq_ss_list, 3)
    assert dict_file.exists()

    kmer_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(kmer_dict[record.id], cgl.get_kmer(record.seq, 3))
            for record in seq_ss_list
        ]
    )


def test_gen_kmer_dict_from_folder(example_id_filter, tmp_path):
    dict_file = tmp_path / "kmer.dict"
    cgl.gen_kmer_dict_from_folder(
        dict_file, cgl.paths.MIRNA_SEQ_SS_FOLDER, 3, keep=example_id_filter
    )
    assert dict_file.exists()

    kmer_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                kmer_dict[id_],
                cgl.get_kmer_by_name(cgl.paths.MIRNA_SEQ_SS_FOLDER, id_, 3),
            )
            for id_ in example_id_filter
        ]
    )
