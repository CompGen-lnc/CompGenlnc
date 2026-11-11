from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from compgenlnc.consts.fields import (
    IS_POSITIVE_FIELD,
    LNCRNA_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
)
from compgenlnc.typing.numpy_dtypes import (
    InteractionPair,
    InteractionTuple,
    dtype_interaction_tuple,
)

if TYPE_CHECKING:
    from _typeshed import StrPath


def save_interactions(
    new_csv: StrPath,
    interactions_list: NDArray[InteractionTuple] | None = None,
    *,
    pos_list: NDArray[InteractionPair] | None = None,
    neg_list: NDArray[InteractionPair] | None = None,
) -> None:
    """Save interactions to a CSV file.

    Parameters
    ----------
    new_csv : str or os.PathLike
        Path to the new CSV file.
    interactions_list : numpy.ndarray of InteractionTuple, optional
        A 1D array with the interaction tuples to save. Cannot be `None`
        if `pos_list` and `neg_list` are not provided. Default: `None`.
    pos_list : numpy.ndarray of InteractionPair, optional
        A 1D array with the positive interaction pairs to save. Default:
        `None`.
    neg_list : numpy.ndarray of InteractionPair, optional
        A 1D array with the negative interaction pairs to save. Default:
        `None`.
    """

    def parse_interactions(interactions, positive):
        """Parse interaction tuples into a structured array."""
        arr = np.empty(interactions.shape, dtype=dtype_interaction_tuple)
        arr[PAIR_FIELDS] = interactions[PAIR_FIELDS]
        arr[IS_POSITIVE_FIELD] = positive
        return arr

    if interactions_list is None:
        pos_interactions = parse_interactions(pos_list, True)
        neg_interactions = parse_interactions(neg_list, False)
        interactions_list = np.concatenate(
            (pos_interactions, neg_interactions)
        )
    df = pd.DataFrame.from_records(interactions_list).sort_values(
        ["positive", "lncRNA", "miRNA"], ascending=[False, True, True]
    )
    df.to_csv(new_csv, index=False)


def load_interactions(
    csv_file: StrPath,
    *,
    lnc_keep: Sequence[str] | None = None,
    mir_keep: Sequence[str] | None = None,
) -> NDArray[InteractionTuple]:
    """Load interactions from a CSV file.

    Parameters
    ----------
    csv_file : str or os.PathLike
        Path of the CSV file with the interactions.
    lnc_keep : sequence of str, optional
        Sequence of lncRNA IDs to filter the interactions. Does not
        apply if `None`.
    mir_keep : sequence of str, optional
        Sequence of miRNA IDs to filter the interactions. Does not
        apply if `None`.

    Returns
    -------
    numpy.ndarray of InteractionTuple
        A 1D array of interaction tuples.
    """
    df = pd.read_csv(csv_file)
    interactions = df.to_records(index=False)
    mask = np.ones(len(interactions), np.bool)
    if lnc_keep is not None and len(lnc_keep):
        mask &= np.isin(interactions[LNCRNA_FIELD], lnc_keep)
    if mir_keep is not None and len(mir_keep):
        mask &= np.isin(interactions[MIRNA_FIELD], mir_keep)
    return interactions[mask]


def load_precursors(precursors_file: StrPath) -> dict[str, str]:
    """Load the mature-precursor miRNA associations.

    Parameters
    ----------
    precursors_file : str or os.PathLike
        Path of the CSV file with the miRNA associations.

    Returns
    -------
    dict of str to str
        A dictionary mapping mature miRNA to its precursor.
    """
    df = pd.read_csv(precursors_file)
    return {mat: pre for mat, pre in df.itertuples(index=False)}
