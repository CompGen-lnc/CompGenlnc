import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            "ACGTGAAC",
            np.array(
                [
                    0.375,
                    0.125,
                    0.25,
                    0.25,
                    0.0,
                    1 / 7,
                    2 / 7,
                    2 / 7,
                    0.0,
                    1 / 7,
                    0.125,
                    0.125,
                    0.75,
                    0.75,
                    0.875,
                    0.5,
                    0.5,
                    0.5,
                    0.5,
                    0.5,
                    0.375,
                    0.375,
                    0.375,
                    0.625,
                    0.625,
                    0.25,
                    0.25,
                    0.25,
                    1,
                    1,
                ],
                np.float32,
            ),
        ),
        (
            "AAA",
            np.array(
                [1] + [0.0] * 9 + [1 / 3, 1 / 3, 2 / 3, 2 / 3, 1] + [0.0] * 15,
                np.float32,
            ),
        ),
        (
            cgl.structs.SeqSSRecord("", "ATCGATCGATCG", ""),
            np.array(
                [0.25] * 4
                + [
                    3 / 11,
                    2 / 11,
                    0.0,
                    0.0,
                    3 / 11,
                    3 / 11,
                    1 / 12,
                    1 / 12,
                    5 / 12,
                    5 / 12,
                    3 / 4,
                    1 / 6,
                    1 / 6,
                    1 / 2,
                    1 / 2,
                    5 / 6,
                    1 / 3,
                    1 / 3,
                    2 / 3,
                    2 / 3,
                    1,
                    1 / 4,
                    1 / 4,
                    7 / 12,
                    7 / 12,
                    11 / 12,
                ],
                np.float32,
            ),
        ),
    ],
)
def test_get_ctd(seq, expected):
    ctd = cgl.get_ctd(seq)
    assert np.array_equal(ctd, expected)


def test_get_ctd_by_name(example_id_list):
    for id_ in example_id_list:
        ctd = cgl.get_ctd_by_name(cgl.paths.MIRNA_SEQ_SS_DIR, id_)
        seq = cgl.read_sequence(cgl.paths.MIRNA_SEQ_SS_DIR, id_)
        expected = cgl.get_ctd(seq)
        assert np.array_equal(ctd, expected)


def test_gen_ctd_dict(seq_ss_list, tmp_path):
    dict_file = tmp_path / "ctd.dict"
    cgl.gen_ctd_dict(dict_file, seq_ss_list)
    assert dict_file.exists()

    ctd_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(ctd_dict[record.id], cgl.get_ctd(record.seq))
            for record in seq_ss_list
        ]
    )


def test_gen_ctd_dict_from_dir(example_id_filter, tmp_path):
    dict_file = tmp_path / "ctd.dict"
    cgl.gen_ctd_dict_from_dir(
        dict_file, cgl.paths.MIRNA_SEQ_SS_DIR, keep=example_id_filter
    )
    assert dict_file.exists()

    ctd_dict = cgl.load_dict(dict_file)
    assert all(
        [
            np.array_equal(
                ctd_dict[id_],
                cgl.get_ctd_by_name(cgl.paths.MIRNA_SEQ_SS_DIR, id_),
            )
            for id_ in example_id_filter
        ]
    )
