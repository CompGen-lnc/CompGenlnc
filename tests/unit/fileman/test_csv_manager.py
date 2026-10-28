import numpy as np
import pandas as pd

from compgenlnc.consts import (
    IS_POSITIVE_FIELD,
    LNCRNA_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
)
from compgenlnc.fileman import (
    load_interactions,
    load_precursors,
    save_interactions,
)

from constants import MIRNA_FOLDER


def test_save_interactions(interaction_set, tmp_path):
    path = tmp_path / "interactions.csv"
    save_interactions(path, interaction_set)
    assert path.exists()
    df = pd.read_csv(path)
    columns = df.columns.to_numpy()
    assert np.isin(
        [LNCRNA_FIELD, MIRNA_FIELD, IS_POSITIVE_FIELD], columns
    ).all()


def test_load_interactions(interaction_set, tmp_path):
    path = tmp_path / "interactions.csv"
    save_interactions(path, interaction_set)
    expected = np.sort(interaction_set, order=PAIR_FIELDS)
    expected = expected[np.argsort(~expected[IS_POSITIVE_FIELD])]
    interactions = load_interactions(path)
    assert np.array_equal(interactions, expected)


def test_load_precursors():
    path = MIRNA_FOLDER / "data_miRNA.csv"
    precursor_dict = load_precursors(path)
    lines = path.read_text().split("\n")
    headers = lines.pop(0)
    assert "Mature" in headers and "Precursor" in headers
    for line, item in zip(lines, precursor_dict.items()):
        expected = tuple(line.split(","))
        assert item == expected
