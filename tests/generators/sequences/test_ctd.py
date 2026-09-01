import numpy as np
import pytest

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.generators.sequences.ctd import (
    gen_ctd_dict,
    gen_ctd_dict_from_folder,
    get_ctd,
    get_ctd_by_name,
)
from compgenlnc.structs import SeqSSRecord
from compgenlnc.utils.dict_manager import load_dict
from compgenlnc.utils.fasta_manager import read_sequence

from constants import SEQ_SS_FOLDER

@pytest.mark.parametrize(
    'seq, expected',
    [
        ('ACGTGAAC', np.array([
            0.375, 0.125, 0.25, 0.25,
            0.0, 1 / 7, 2 / 7, 2 / 7, 0.0, 1 / 7,
            0.125, 0.125, 0.75, 0.75, 0.875, 
            0.5, 0.5, 0.5, 0.5, 0.5,
            0.375, 0.375, 0.375, 0.625, 0.625,
            0.25, 0.25, 0.25, 1, 1,
        ], np.float32)),
        ('AAA', np.array(
            [1] + [0.0] * 9
            + [1 / 3, 1 / 3, 2 / 3, 2 / 3, 1]
            + [0.0] * 15,
            np.float32
        )),
        (
            SeqSSRecord('', 'ATCGATCGATCG', ''),
            np.array([0.25]*4 + [
                3 / 11, 2 / 11, 0.0, 0.0, 3 / 11, 3 / 11,
                1 / 12, 1 / 12, 5 / 12, 5 / 12, 3 / 4,
                1 / 6, 1 / 6, 1 / 2, 1 / 2, 5 / 6,
                1 / 3, 1 / 3, 2 / 3, 2 / 3, 1,
                1 / 4, 1 / 4, 7 / 12, 7 / 12, 11 / 12,
            ], np.float32)
        )
    ]
)
def test_get_ctd(seq, expected):
    ctd = get_ctd(seq)
    assert np.array_equal(ctd, expected)


def test_get_ctd_by_name(mirna_id):
    ctd = get_ctd_by_name(mirna_id, SEQ_SS_FOLDER)
    seq = read_sequence(SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{mirna_id}.dat')
    expected = get_ctd(seq)
    assert np.array_equal(ctd, expected)


def test_gen_ctd_dict(seq_ss_filtered_list, tmp_path):
    dict_file = tmp_path / 'ctd.dict'
    gen_ctd_dict(seq_ss_filtered_list, dict_file)
    assert dict_file.exists()

    ctd_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            ctd_dict[record.id],
            get_ctd(record.seq)
        )
        for record in seq_ss_filtered_list
    ])


def test_gen_ctd_dict_from_folder(mirna_id_filter, tmp_path):
    dict_file = tmp_path / 'ctd.dict'
    gen_ctd_dict_from_folder(SEQ_SS_FOLDER, dict_file, keep=mirna_id_filter)
    assert dict_file.exists()

    ctd_dict = load_dict(dict_file)
    assert all([
        np.array_equal(
            ctd_dict[id_],
            get_ctd_by_name(id_, SEQ_SS_FOLDER)
        )
        for id_ in mirna_id_filter
    ])