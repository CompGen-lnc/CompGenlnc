import os

import numpy as np
import pandas as pd

from compgenlnc.typing.numpy_dtypes import (
    dtype_interaction_tuple,
    InteractionPair,
    InteractionTuple,
)


def save_interactions(
    new_csv: str | os.PathLike,
    interactions_list: np.typing.NDArray[InteractionTuple] | None = None,
    *,
    pos_list: np.typing.NDArray[InteractionPair] | None = None,
    neg_list: np.typing.NDArray[InteractionPair] | None = None,
) -> None:
    parse_interactions = lambda interactions, positive: np.full(
        (interactions.shape[0], interactions.shape[1] + 1),
        positive,
        dtype_interaction_tuple,
    )

    if interactions_list is None:
        pos_interactions = parse_interactions(pos_list, True)
        neg_interactions = parse_interactions(neg_list, False)
        interactions = np.stack((pos_interactions, neg_interactions))
    df = pd.DataFrame.from_records(interactions).sort_values(
        ["positive", "lncRNA", "miRNA"], ascending=[False, True, True]
    )
    df.to_csv(new_csv, index=False)


def load_interactions(
    filename: str | os.PathLike,
) -> np.typing.NDArray[InteractionTuple]:
    df = pd.read_csv(filename)
    return np.fromiter(
        df.itertuples(index=False), dtype_interaction_tuple, len(df.index)
    )


def load_precursors(precursors_file: str | os.PathLike) -> dict[str, str]:
    df = pd.read_csv(precursors_file)
    return {mat: pre for mat, pre in df.itertuples(index=False)}
