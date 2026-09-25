import os

import numpy as np
import pandas as pd

from compgenlnc.typing.numpy_dtypes import (
    dtype_interaction_tuple,
    InteractionPair,
    InteractionTuple,
)


def save_interactions(
    new_csv: os.PathLike,
    interactions_list: np.typing.NDArray[InteractionTuple] | None = None,
    *,
    pos_list: np.typing.NDArray[InteractionPair] | None = None,
    neg_list: np.typing.NDArray[InteractionPair] | None = None,
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


def load_interactions(
    filename: os.PathLike,
) -> np.typing.NDArray[InteractionTuple]:
    df = pd.read_csv(filename)
    return np.fromiter(
        df.itertuples(index=False), dtype_interaction_tuple, len(df.index)
    )
