import os

import numpy as np

import compgenlnc as cgl


def test_get_mature_index(interaction_list, precursors_consverter, subtests):
    for interaction in interaction_list[cgl.PAIR_FIELDS]:
        index = cgl.get_mature_index(interaction, precursors_consverter)
        lnc, mir = interaction
        with subtests.test(i=(lnc, mir)):
            if mir not in precursors_consverter:
                assert index is None
            else:
                assert index[cgl.LNCRNA_FIELD] == lnc
                assert index[cgl.MATURE_FIELD] == mir
                assert index[cgl.PRECURSOR_FIELD] == precursors_consverter[mir]


def test_gen_indexes(interaction_set, precursors_consverter, tmp_path):
    path = tmp_path / "index.npy"
    cgl.gen_indexes(interaction_set, precursors_consverter, index_file=path)
    assert path.exists()
    indexes = np.load(path, allow_pickle=True)
    for lnc, prec, mat, inter in indexes:
        mask: np.ndarray = (interaction_set[cgl.LNCRNA_FIELD] == lnc) & (
            interaction_set[cgl.MIRNA_FIELD] == mat
        )
        assert mask.any()
        pair = interaction_set[mask][0]
        assert prec == precursors_consverter[pair[cgl.MIRNA_FIELD]]
        assert inter == f"{lnc}_{prec}" or inter == f"{lnc}_{mat}"


def test_load_indexes():
    path = cgl.paths.EMBEDDING_DIR / "index.npy"
    indexes = cgl.load_indexes(index_file=path)
    array = np.load(path, allow_pickle=True)
    expected = np.array(array, dtype=cgl.dtype_index)
    assert np.array_equal(indexes, expected)


def test_dict_to_numpy(index_list, tmp_path, subtests):
    dir = cgl.paths.MAT_MIRNA_DICTS_DIR
    index_path = cgl.paths.EMBEDDING_DIR / "index.npy"
    dicts = os.listdir(dir)
    for filename in dicts:
        dict_path = dir / filename
        npy_path = tmp_path / filename.replace(".dict", ".npy")
        with subtests.test(i=(filename)):
            dict_ = cgl.load_dict(dict_path)
            cgl.dict_to_numpy(
                npy_path, dict_, index_file=index_path, rna=cgl.MATURE
            )
            feat_array = np.load(npy_path)
            for index, vector in zip(index_list, feat_array):
                assert np.array_equal(dict_[index[cgl.MATURE_FIELD]], vector)
