import os

import numpy as np

from compgenlnc import (
    dict_to_numpy,
    gen_indexes,
    get_mature_index,
    load_indexes,
)
from compgenlnc.consts import (
    LNCRNA_FIELD,
    MATURE,
    MATURE_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.fileman import load_dict
from compgenlnc.typing import dtype_index

from constants import EMBEDDING_FOLDER, MIRNA_FOLDER


def test_get_mature_index(interaction_list, precursors_consverter, subtests):
    for interaction in interaction_list[PAIR_FIELDS]:
        index = get_mature_index(interaction, precursors_consverter)
        lnc, mir = interaction
        with subtests.test(i=(lnc, mir)):
            if mir not in precursors_consverter:
                assert index == None
            else:
                assert index[LNCRNA_FIELD] == lnc
                assert index[MATURE_FIELD] == mir
                assert index[PRECURSOR_FIELD] == precursors_consverter[mir]


def test_gen_indexes(interaction_set, precursors_consverter, tmp_path):
    path = tmp_path / "index.npy"
    gen_indexes(path, interaction_set, precursors_consverter)
    assert path.exists()
    indexes = np.load(path, allow_pickle=True)
    for lnc, prec, mat in indexes:
        mask: np.ndarray = (interaction_set[LNCRNA_FIELD] == lnc) & (
            interaction_set[MIRNA_FIELD] == mat
        )
        assert mask.any()
        pair = interaction_set[mask][0]
        assert prec == precursors_consverter[pair[MIRNA_FIELD]]


def test_load_indexes():
    path = EMBEDDING_FOLDER / "index.npy"
    indexes = load_indexes(path)
    array = np.load(path, allow_pickle=True)
    expected = np.array(array, dtype=dtype_index)
    assert np.array_equal(indexes, expected)


def test_dict_to_numpy(index_list, tmp_path, subtests):
    folder = MIRNA_FOLDER / "dict" / "mature"
    dicts = os.listdir(folder)
    for filename in dicts:
        dict_path = folder / filename
        npy_path = tmp_path / filename.replace(".dict", ".npy")
        with subtests.test(i=(filename)):
            dict_ = load_dict(dict_path)
            dict_to_numpy(npy_path, dict_, index_list, rna=MATURE)
            feat_array = np.load(npy_path, allow_pickle=True)
            for index, vector in zip(index_list, feat_array):
                assert np.array_equal(dict_[index[MATURE_FIELD]], vector)