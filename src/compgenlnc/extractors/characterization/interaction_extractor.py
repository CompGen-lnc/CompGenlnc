import os
import random
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd


dtype_interaction_tuple = np.dtype(
    [
        ("lncRNA", np.dtypes.StringDType),
        ("miRNA", np.dtypes.StringDType),
        ("positive", np.bool),
    ]
)

dtype_interaction_pair = np.dtype(
    [
        ("lncRNA", np.dtypes.StringDType),
        ("miRNA", np.dtypes.StringDType),
    ]
)


def save_interactions(
    new_csv: os.PathLike,
    interactions_list: np.typing.NDArray[np.void] | None = None,
    *,
    pos_list: np.typing.NDArray[np.void] | None = None,
    neg_list: np.typing.NDArray[np.void] | None = None,
) -> None:
    parse_interactions = lambda interactions, positive: {
        "lncRNA": interactions["lncRNA"],
        "miRNA": interactions["miRNA"],
        "positive": np.array([positive] * len(interactions), dtype=np.bool),
    }

    if interactions_list is not None:
        interactions = {
            "lncRNA": interactions_list["lncRNA"],
            "miRNA": interactions_list["miRNA"],
            "positive": interactions_list["positive"],
        }
    else:
        pos_interactions = parse_interactions(pos_list, True)
        neg_interactions = parse_interactions(neg_list, False)
        interactions = {
            key: np.concatenate((pos_interactions[key], neg_interactions[key]))
            for key in pos_interactions.keys()
        }
    df = pd.DataFrame(interactions).sort_values(
        ["positive", "lncRNA", "miRNA"], ascending=[False, True, True]
    )
    df.to_csv(new_csv, index=False)


def load_interactions(filename: os.PathLike) -> np.typing.NDArray[np.void]:
    df = pd.read_csv(filename)
    return np.fromiter(
        df.itertuples(index=False), dtype_interaction_tuple, len(df.index)
    )


def extract_interactions(
    new_csv: os.PathLike,
    pos_file: os.PathLike,
    neg_file: os.PathLike | None = None,
    /,
    separator: str | None = None,
    lnc_list: Sequence[str] | None = None,
    mir_list: Sequence[str] | None = None,
) -> None:
    if neg_file is None and (not lnc_list or not mir_list):
        raise ValueError(
            "If neg_file is None, "
            "lnc_list nor mir_list cannot be None or empty"
        )

    def read_interactions(file_inter: Path) -> np.typing.NDArray[np.void]:
        interactions_list = []
        if not file_inter:
            return interactions_list
        with open(file_inter) as file:
            for line in file.readlines():
                line = line.strip()
                if not line or line == "lncRNA,miRNA":
                    continue
                interaction = line.split(separator)
                interactions_list.append(
                    tuple((id_.strip() for id_ in interaction))
                )

        print(interactions_list)
        return np.array(interactions_list, dtype_interaction_pair)

    pos_list = read_interactions(pos_file)
    if neg_file is not None:
        neg_list = read_interactions(neg_file)
    else:
        length = len(pos_list)
        lnc_neg_list = random.choices(lnc_list, k=length)
        mir_neg_list = random.choices(mir_list, k=length)
        for i in range(length):
            while True:
                pair = (pos_list["lncRNA"] == lnc_neg_list[i]) & (
                    pos_list["miRNA"] == mir_neg_list[i]
                )
                if not pair.any():
                    break
                lnc_neg_list[i] = random.choice(lnc_list)
                mir_neg_list[i] = random.choice(mir_list)
        neg_interactions = set(zip(lnc_neg_list, mir_neg_list))
        neg_list = np.array(list(neg_interactions), dtype_interaction_pair)
    save_interactions(new_csv, pos_list=pos_list, neg_list=neg_list)


def filter_interactions(
    new_csv: os.PathLike,
    interactions_file: os.PathLike,
    /,
    lnc_filter: Sequence[str] | None = None,
    mir_filter: Sequence[str] | None = None,
) -> None:
    interactions = load_interactions(interactions_file)
    if lnc_filter:
        interactions = interactions[
            np.isin(interactions["lncRNA"], lnc_filter)
        ]
    if mir_filter:
        interactions = interactions[np.isin(interactions["miRNA"], mir_filter)]
    save_interactions(new_csv, interactions)
