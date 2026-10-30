import numpy as np
import pandas as pd

import compgenlnc as cgl


def test_save_interactions(interaction_set, tmp_path):
    path = tmp_path / "interactions.csv"
    cgl.save_interactions(path, interaction_set)
    assert path.exists()
    df = pd.read_csv(path)
    columns = df.columns.to_numpy()
    assert np.isin(
        [cgl.LNCRNA_FIELD, cgl.MIRNA_FIELD, cgl.IS_POSITIVE_FIELD], columns
    ).all()


def test_load_interactions(interaction_set, tmp_path):
    path = tmp_path / "interactions.csv"
    cgl.save_interactions(path, interaction_set)
    positives = interaction_set[interaction_set[cgl.IS_POSITIVE_FIELD]]
    negatives = interaction_set[~interaction_set[cgl.IS_POSITIVE_FIELD]]
    positives.sort(order=cgl.PAIR_FIELDS)
    negatives.sort(order=cgl.PAIR_FIELDS)
    expected = np.concatenate((positives, negatives), axis=None)
    interactions = cgl.load_interactions(path)
    assert np.array_equal(interactions, expected)


def test_load_precursors():
    path = cgl.paths.MIRNA_FOLDER / "data_miRNA.csv"
    precursor_dict = cgl.load_precursors(path)
    lines = path.read_text().split("\n")
    headers = lines.pop(0)
    assert "Mature" in headers and "Precursor" in headers
    for line, item in zip(lines, precursor_dict.items()):
        expected = tuple(line.split(","))
        assert item == expected
