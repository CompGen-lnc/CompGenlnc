import os
from collections.abc import Iterable

import numpy as np

from compgenlnc.consts.fields import (
    LNCRNA_FIELD,
    MATURE_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.consts.params import LNCRNA, MATURE, PRECURSOR
from compgenlnc.typing.numpy_dtypes import dtype_index, Index, InteractionPair


def get_mature_index(
    interaction: InteractionPair, precursor_dict: dict[str, str]
) -> Index:
    lnc, mir = interaction[PAIR_FIELDS]
    if mir not in precursor_dict.keys():
        return None
    precursor = precursor_dict[mir]
    return np.array((lnc, precursor, mir), dtype=dtype_index)


def gen_indexes(
    index_file: str | os.PathLike,
    interactions_it: Iterable[InteractionPair],
    precursor_dict: dict[str, str],
) -> None:
    interactions_it = list(interactions_it)
    matures = precursor_dict.keys()
    array = np.empty((len(interactions_it), 3), dtype=dtype_index)
    fails = []
    for i, interaction in enumerate(interactions_it):
        mir = interaction[MIRNA_FIELD]
        if mir in matures:
            array[i] = get_mature_index(interaction, precursor_dict)
        else:
            fails.append(i)
    mask = np.ones(len(array), dtype=bool)
    mask[fails] = False
    indexes = np.unique(array[mask])
    np.save(index_file, indexes)


def load_indexes(index_file: str | os.PathLike) -> np.typing.NDArray[Index]:
    return np.load(index_file, allow_pickle=True)


def dict_to_numpy(
    np_file: str | os.PathLike,
    feat_dict: dict[str, np.typing.NDArray],
    indexes: np.typing.NDArray[Index],
    *,
    rna: str = LNCRNA,
) -> None:
    if rna == LNCRNA:
        index_rna = indexes[LNCRNA_FIELD]
    elif rna == MATURE:
        index_rna = indexes[MATURE_FIELD]
    elif rna == PRECURSOR:
        index_rna = indexes[PRECURSOR_FIELD]

    sample_vec = next(iter(feat_dict.values()))
    array = np.empty((len(index_rna), len(sample_vec)), dtype=sample_vec.dtype)
    for i, index in enumerate(index_rna):
        if index in feat_dict.keys():
            array[i] = feat_dict[index]
        else:
            array[i] = 0
    np.save(np_file, array)
