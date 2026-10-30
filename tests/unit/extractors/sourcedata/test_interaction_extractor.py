import numpy as np

import compgenlnc as cgl


def test_save_interactions(
    interaction_set,
    interaction_positive_pairs,
    interaction_negative_pairs,
    tmp_path,
):
    filename = tmp_path / "interaction.csv"
    cgl.save_interactions(filename, interaction_set)
    assert filename.exists()
    joint_filename = tmp_path / "joint_interaction.csv"
    cgl.save_interactions(
        joint_filename,
        pos_list=interaction_positive_pairs,
        neg_list=interaction_negative_pairs,
    )
    assert joint_filename.exists()
    assert filename.read_text() == joint_filename.read_text()

    lines = filename.read_text().split("\n")
    for line, expected in zip(lines[1:], interaction_set):
        values = line.split(",")
        values[2] = values[2] == "True"
        arr = np.array(tuple(values), cgl.dtype_interaction_tuple)
        assert np.array_equal(arr, expected)


def test_load_interactions(interaction_set, tmp_path):
    filename = tmp_path / "interaction.csv"
    cgl.save_interactions(filename, interaction_set)
    loaded_list = cgl.load_interactions(filename)
    assert np.array_equal(loaded_list, interaction_set)


def test_extract_interactions_without_negavites(
    lncrna_id_list, mirna_id_list, tmp_path
):
    filename = tmp_path / "interaction.csv"
    cgl.extract_interactions(
        filename,
        cgl.paths.PAIRS_FILE_POSITIVE,
        separator=",",
        lnc_list=lncrna_id_list,
        mir_list=mirna_id_list,
    )
    assert filename.exists()


def test_extract_interactions_with_negavites(tmp_path):
    filename = tmp_path / "interaction.csv"
    cgl.extract_interactions(
        filename,
        cgl.paths.PAIRS_FILE_POSITIVE,
        cgl.paths.PAIRS_FILE_NEGATIVE,
        separator=",",
    )
    assert filename.exists()
    assert np.array_equal(
        cgl.load_interactions(filename),
        cgl.load_interactions(cgl.paths.PAIRS_FILE),
    )


def test_filter_interactions(lncrna_id_filter, mirna_id_filter, tmp_path):
    filename = tmp_path / "interaction.csv"
    cgl.filter_interactions(
        filename,
        cgl.paths.PAIRS_FILE,
        lnc_filter=lncrna_id_filter,
        mir_filter=mirna_id_filter,
    )
    arr = cgl.load_interactions(filename)
    expected = cgl.load_interactions(cgl.paths.PAIRS_FILE)
    mask = expected[cgl.IS_POSITIVE_FIELD] == expected[cgl.IS_POSITIVE_FIELD]
    if lncrna_id_filter:
        mask &= np.isin(expected[cgl.LNCRNA_FIELD], lncrna_id_filter)
    if mirna_id_filter:
        mask &= np.isin(expected[cgl.MIRNA_FIELD], mirna_id_filter)
    assert np.array_equal(arr, expected[mask])
