import numpy as np

from compgenlnc.extractors.interactions import (
    extract_interactions,
    filter_interactions,
)
from compgenlnc.typing import dtype_interaction_tuple
from compgenlnc.utils.interactions_manager import (
    load_interactions,
    save_interactions,
)

from constants import (
    INTERACTIONS_CSV,
    NEG_INTERACTIONS_CSV,
    POS_INTERACTIONS_CSV,
)


def test_save_interactions(
    interaction_set,
    interaction_positive_pairs,
    interaction_negative_pairs,
    tmp_path,
):
    filename = tmp_path / "interaction.csv"
    save_interactions(filename, interaction_set)
    assert filename.exists()
    joint_filename = tmp_path / "joint_interaction.csv"
    save_interactions(
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
        arr = np.array(tuple(values), dtype_interaction_tuple)
        assert np.array_equal(arr, expected)


def test_load_interactions(interaction_set, tmp_path):
    filename = tmp_path / "interaction.csv"
    save_interactions(filename, interaction_set)
    loaded_list = load_interactions(filename)
    assert np.array_equal(loaded_list, interaction_set)


def test_extract_interactions_without_negavites(
    lncrna_id_list, mirna_id_list, tmp_path
):
    filename = tmp_path / "interaction.csv"
    extract_interactions(
        filename,
        POS_INTERACTIONS_CSV,
        separator=",",
        lnc_list=lncrna_id_list,
        mir_list=mirna_id_list,
    )
    assert filename.exists()


def test_extract_interactions_with_negavites(tmp_path):
    filename = tmp_path / "interaction.csv"
    extract_interactions(
        filename, POS_INTERACTIONS_CSV, NEG_INTERACTIONS_CSV, separator=","
    )
    assert filename.exists()
    assert np.array_equal(
        load_interactions(filename), load_interactions(INTERACTIONS_CSV)
    )


def test_filter_interactions(lncrna_id_filter, mirna_id_filter, tmp_path):
    filename = tmp_path / "interaction.csv"
    filter_interactions(
        filename,
        INTERACTIONS_CSV,
        lnc_filter=lncrna_id_filter,
        mir_filter=mirna_id_filter,
    )
    arr = load_interactions(filename)
    expected = load_interactions(INTERACTIONS_CSV)
    mask = expected["positive"] == expected["positive"]
    if lncrna_id_filter:
        mask &= np.isin(expected["lncRNA"], lncrna_id_filter)
    if mirna_id_filter:
        mask &= np.isin(expected["miRNA"], mirna_id_filter)
    assert np.array_equal(arr, expected[mask])
