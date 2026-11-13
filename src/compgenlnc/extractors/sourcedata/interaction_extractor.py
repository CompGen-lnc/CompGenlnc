from __future__ import annotations

import random
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from compgenlnc.consts.fields import LNCRNA_FIELD, MIRNA_FIELD
from compgenlnc.dtypes import dtype_interaction_pair
from compgenlnc.fileman.csv_manager import (
    load_interactions,
    save_interactions,
)
from compgenlnc.typing.numpy_dtypes import (
    InteractionPair,
)

if TYPE_CHECKING:
    from _typeshed import StrPath


def extract_interactions(
    new_csv: StrPath,
    pos_file: StrPath,
    neg_file: StrPath | None = None,
    *,
    separator: str | None = None,
    lnc_list: Sequence[str] | None = None,
    mir_list: Sequence[str] | None = None,
) -> None:
    """Extract the interactions from the positive and negative files.

    Parameters
    ----------
    new_csv : str or os.PathLike
        Path of the CSV file where the interactions will be saved.
    pos_file : str or os.PathLike
        Path of the file containing the positive interactions.
    neg_file : str or os.PathLike, optional
        Path of the file containing the negative interactions.
        Default: `None`.
    separator : str, optional
        Separator used in the paositive and negative files. Default:
        `None`.
    lnc_list : sequence of str, optional
        Sequence of lncRNA IDs to be included in the interactions. Does
        not apply if `None` or empty. Default: `None`.
    mir_list : sequence of str, optional
        Sequence of lncRNA IDs to be included in the interactions. Does
        not apply if `None` or empty. Default: `None`.
    """
    if neg_file is None and (not lnc_list or not mir_list):
        raise ValueError(
            "If neg_file is None, "
            "lnc_list nor mir_list cannot be None or empty"
        )

    def read_interactions(
        file_inter: StrPath,
    ) -> NDArray[InteractionPair]:
        """Read interactions from the file."""
        interactions_list = []
        if not file_inter:
            return interactions_list
        with open(file_inter) as file:
            for line in file.readlines():
                line = line.strip()
                if not line or "lncRNA" in line or "miRNA" in line:
                    continue
                lnc, mir = line.split(separator)[:2]
                if (
                    lnc_list is not None
                    and lnc not in lnc_list
                    or mir_list is not None
                    and mir not in mir_list
                ):
                    continue
                interactions_list.append((lnc.strip(), mir.strip()))

        return np.unique(np.array(interactions_list, dtype_interaction_pair))

    pos_list = read_interactions(pos_file)
    if neg_file is not None and Path(neg_file).resolve().exists():
        neg_list = read_interactions(neg_file)
    else:
        length = len(pos_list)
        lnc_neg_list = random.choices(lnc_list, k=length)
        mir_neg_list = random.choices(mir_list, k=length)
        for i in range(length):
            while True:
                pair = (pos_list[LNCRNA_FIELD] == lnc_neg_list[i]) & (
                    pos_list[MIRNA_FIELD] == mir_neg_list[i]
                )
                if not pair.any():
                    break
                lnc_neg_list[i] = random.choice(lnc_list)
                mir_neg_list[i] = random.choice(mir_list)
        neg_interactions = set(zip(lnc_neg_list, mir_neg_list))
        neg_list = np.array(list(neg_interactions), dtype_interaction_pair)
        df = pd.DataFrame.from_records(neg_list)
        if neg_file is not None:
            df.to_csv(neg_file, index=False)
    save_interactions(new_csv, pos_list=pos_list, neg_list=neg_list)


def filter_interactions(
    interactions_file: StrPath,
    new_csv: StrPath,
    /,
    lnc_filter: Sequence[str] | None = None,
    mir_filter: Sequence[str] | None = None,
) -> None:
    """Filter interactions based on the provided filters.

    Parameters
    ----------
    interactions_file : str or os.PathLike
        Path of the CSV file with interactions.
    new_csv : str or os.PathLike
        Path of the CSV file to save filtered interactions.
    lnc_filter : sequence of str, optional
        List of lncRNAs to filter interactions by.
    mir_filter : sequence of str, optional
        List of miRNAs to filter interactions by.
    """
    interactions = load_interactions(interactions_file)
    if lnc_filter:
        interactions = interactions[
            np.isin(interactions[LNCRNA_FIELD], lnc_filter)
        ]
    if mir_filter:
        interactions = interactions[
            np.isin(interactions[MIRNA_FIELD], mir_filter)
        ]
    save_interactions(new_csv, interactions)
