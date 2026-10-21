import os
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from compgenlnc.consts.paths import (
    INDEX_NPY,
    LNCRNA_SEQ_SS_FOLDER,
    MAT_MIRNA_FASTA,
    MIRNA_SEQ_SS_FOLDER,
)
from compgenlnc.consts.fields import (
    INTERACTION_FIELD,
    LNCRNA_FIELD,
    MATURE_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.consts.params import INTERACTION, LNCRNA, MATURE, PRECURSOR
from compgenlnc.fileman.csv_manager import load_interactions, load_precursors
from compgenlnc.fileman.fasta_manager import iter_seq_ss, load_fasta
from compgenlnc.typing.numpy_dtypes import dtype_index, Index, InteractionPair


_cgl_indexes = None
_lnc_indexed_records = None
_pre_mir_indexed_records = None
_mat_mir_indexed_records = None


def get_mature_index(
    interaction: InteractionPair, precursor_dict: dict[str, str]
) -> Index:
    lnc, mir = interaction[PAIR_FIELDS]
    if mir not in precursor_dict.keys():
        return None
    precursor = precursor_dict[mir]
    inter_str = f"{lnc}_{mir}"
    return np.array((lnc, precursor, mir, inter_str), dtype=dtype_index)


def gen_indexes(
    interactions_it: Iterable[InteractionPair],
    precursor_dict: dict[str, str],
    *,
    index_file: str | os.PathLike = INDEX_NPY,
) -> None:
    interactions_it = list(interactions_it)
    array = np.empty((len(interactions_it), 3), dtype=dtype_index)
    fails = []
    for i, interaction in enumerate(interactions_it):
        mir = interaction[MIRNA_FIELD]
        if mir in precursor_dict:
            array[i] = get_mature_index(interaction, precursor_dict)
        else:
            fails.append(i)
    mask = np.ones(len(array), dtype=bool)
    mask[fails] = False
    masked_array = array[mask]
    unique_array, order = np.unique(masked_array, return_index=True)
    _cgl_indexes = unique_array[np.argsort(order)]
    index_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(index_file, _cgl_indexes)


def gen_indexes_from_files(
    interaction_csv: str | os.PathLike,
    precursor_csv: str | os.PathLike,
    *,
    lnc_list: Iterable[str] | None = None,
    mir_list: Iterable[str] | None = None,
) -> None:
    interactions = load_interactions(interaction_csv)
    mask = np.ones(len(interactions), dtype=np.bool)
    if lnc_list is not None and len(lnc_list):
        mask &= np.isin(interactions[LNCRNA_FIELD], list(lnc_list))
    if mir_list is not None and len(mir_list):
        mask &= np.isin(interactions[MIRNA_FIELD], list(mir_list))
    interactions = interactions[mask]
    precursors = load_precursors(precursor_csv)
    gen_indexes(interactions, precursors)


def load_indexes(
    *, index_file: str | os.PathLike = INDEX_NPY, rna: str | None = None
) -> np.typing.NDArray[Index]:
    global _cgl_indexes
    if _cgl_indexes is None:
        _cgl_indexes = np.load(index_file, allow_pickle=True)

    if rna == LNCRNA:
        return _cgl_indexes[LNCRNA_FIELD]
    elif rna == MATURE:
        return _cgl_indexes[MATURE_FIELD]
    elif rna == PRECURSOR:
        return _cgl_indexes[PRECURSOR_FIELD]
    elif rna == INTERACTION:
        return _cgl_indexes[INTERACTION_FIELD]
    return _cgl_indexes


def index_lncrna(
    *, index_file: str | os.PathLike = INDEX_NPY
) -> np.typing.NDArray[np.object_]:
    global _lnc_indexed_records
    if _lnc_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=LNCRNA)
        records = {
            rec.id: rec
            for rec in iter_seq_ss(LNCRNA_SEQ_SS_FOLDER, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _lnc_indexed_records = map_index(indexes)
    return _lnc_indexed_records


def index_precursor_mirna(
    *, index_file: str | os.PathLike = INDEX_NPY
) -> np.typing.NDArray[np.object_]:
    global _pre_mir_indexed_records
    if _pre_mir_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=PRECURSOR)
        records = {
            rec.id: rec
            for rec in iter_seq_ss(MIRNA_SEQ_SS_FOLDER, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _pre_mir_indexed_records = map_index(indexes)
    return _pre_mir_indexed_records


def index_mature_mirna(
    *, index_file: str | os.PathLike = INDEX_NPY
) -> np.typing.NDArray[np.object_]:
    global _mat_mir_indexed_records
    if _mat_mir_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=MATURE)
        records = {
            rec.id: rec for rec in load_fasta(MAT_MIRNA_FASTA, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _mat_mir_indexed_records = map_index(indexes)
    return _mat_mir_indexed_records


def dict_to_numpy(
    np_file: str | os.PathLike,
    feat_dict: dict[str, np.typing.NDArray],
    *,
    index_file: str | os.PathLike = INDEX_NPY,
    rna: str = LNCRNA,
) -> None:
    indexes = load_indexes(index_file=index_file, rna=rna)

    sample_vec = next(iter(feat_dict.values()))
    array = np.empty((len(indexes), len(sample_vec)), dtype=sample_vec.dtype)
    for i, index in enumerate(indexes):
        if index in feat_dict.keys():
            array[i] = feat_dict[index]
        else:
            array[i] = 0
    np_file = Path(np_file).resolve()
    np_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(np_file, array)
