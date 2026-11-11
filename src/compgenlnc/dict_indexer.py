from __future__ import annotations

from collections.abc import Collection, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from compgenlnc.consts.fields import (
    INTERACTION_FIELD,
    LNCRNA_FIELD,
    MATURE_FIELD,
    MIRNA_FIELD,
    PAIR_FIELDS,
    PRECURSOR_FIELD,
)
from compgenlnc.consts.params import INTERACTION, LNCRNA, MATURE, PRECURSOR
from compgenlnc.consts.paths import (
    INDEX_NPY,
    LNCRNA_SEQ_SS_DIR,
    MAT_MIRNA_FASTA,
    MIRNA_SEQ_SS_DIR,
    PAIRS_FILE,
)
from compgenlnc.fileman.csv_manager import load_interactions, load_precursors
from compgenlnc.fileman.fasta_manager import iter_fasta, iter_seq_ss
from compgenlnc.typing.numpy_dtypes import (
    Index,
    InteractionPair,
    dtype_index,
    dtype_interaction_pair,
)

if TYPE_CHECKING:
    from _typeshed import StrPath

_cgl_indexes = None
_index_file = None
_lnc_indexed_records = None
_pre_mir_indexed_records = None
_mat_mir_indexed_records = None
_indexed_interactions = None


def get_mature_index(
    interaction: InteractionPair, precursor_dict: dict[str, str]
) -> Index:
    """Get the index for an lncRNA-mature miRNA interaction.

    Parameters
    ----------
    interaction : InteractionPair
        The interaction to get the index for.
    precursor_dict : dict of str to str
        Mapping of the mature miRNA to its precursor.

    Returns
    -------
    Index
        The index for the interaction. If the miRNA is not in
        `precursor_dict`, `None` is returned.
    """
    lnc, mir = interaction[PAIR_FIELDS]
    if mir not in precursor_dict.keys():
        return None
    precursor = precursor_dict[mir]
    inter_str = f"{lnc}_{mir}"
    return np.array((lnc, precursor, mir, inter_str), dtype=dtype_index)


def gen_indexes(
    interactions: Collection[InteractionPair],
    precursor_dict: dict[str, str],
    *,
    index_file: StrPath = INDEX_NPY,
) -> None:
    """Generate the indexes from an interaction set.

    Parameters
    ----------
    interactions : collection of InteractionPair
        Set of interaction pairs to index.
    precursor_dict : dict of str to str
        Mapping of the mature miRNA to its precursor.
    index_file : str or os.PathLike, optional
        Path of the file to save the indexes. Default:
        `paths.INDEX_NPY`.
    """
    array = np.empty((len(interactions), 3), dtype=dtype_index)
    fails = []
    for i, interaction in enumerate(interactions):
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
    Path(index_file).parent.mkdir(parents=True, exist_ok=True)
    np.save(index_file, _cgl_indexes)


def gen_indexes_from_files(
    interaction_csv: StrPath,
    precursor_csv: StrPath,
    *,
    lnc_list: Sequence[str] | None = None,
    mir_list: Sequence[str] | None = None,
) -> None:
    """Generate the indexes from interaction and precursor files.

    Parameters
    ----------
    interaction_csv : str or os.PathLike
        Path of the CSV file with the interaction pairs.
    precursor_csv : str or os.PathLike
        Path of the CSV file with the precursor-mature miRNA pairs.
    lnc_list : sequence of str, optional
        Set of IDs of LNCRNAs to include in the indexes.
    mir_list : sequence of str, optional
        Set of IDs of miRNAs to include in the indexes.
    """
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
    *, index_file: StrPath = INDEX_NPY, rna: str | None = None
) -> NDArray[Index]:
    """Load the indexes from a file.

    Parameters
    ----------
    index_file : str or os.PathLike, optional
        Path of the file containing the indexes. Default:
        `paths.INDEX_PY`.
    rna : str, optional
        The index column to return. If `None`, all columns are returned.
        Default: `None`.

    Returns
    -------
    numpy.ndarray of Index
        A 1D array with indexes for the database.
    """
    global _cgl_indexes
    global _lnc_indexed_records
    global _pre_mir_indexed_records
    global _mat_mir_indexed_records
    global _indexed_interactions
    global _index_file
    if _cgl_indexes is None or _index_file != index_file:
        _cgl_indexes = np.load(index_file, allow_pickle=True)
        _lnc_indexed_records = None
        _pre_mir_indexed_records = None
        _mat_mir_indexed_records = None
        _indexed_interactions = None
        _index_file = index_file

    if rna == LNCRNA:
        return _cgl_indexes[LNCRNA_FIELD]
    elif rna == MATURE:
        return _cgl_indexes[MATURE_FIELD]
    elif rna == PRECURSOR:
        return _cgl_indexes[PRECURSOR_FIELD]
    elif rna == INTERACTION:
        return _cgl_indexes[INTERACTION_FIELD]
    return _cgl_indexes


def index_lncrna(*, index_file: StrPath = INDEX_NPY) -> NDArray[np.object_]:
    """Index the lncRNA records.

    Parameters
    ----------
    index_file : str or os.PathLike, optional
        Path of the file containing the indexes. Default:
        `paths.INDEX_PY`.

    Returns
    -------
    numpy.ndarray of SeqSSRecord
        A 1D array of the lncRNA records sorted by indexes.
    """
    global _lnc_indexed_records
    if _lnc_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=LNCRNA)
        records = {
            rec.id: rec for rec in iter_seq_ss(LNCRNA_SEQ_SS_DIR, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _lnc_indexed_records = map_index(indexes)
    return _lnc_indexed_records


def index_precursor_mirna(
    *, index_file: StrPath = INDEX_NPY
) -> NDArray[np.object_]:
    """Index the precursor miRNA records.

    Parameters
    ----------
    index_file : str or os.PathLike, optional
        Path of the file containing the indexes. Default:
        `paths.INDEX_PY`.

    Returns
    -------
    numpy.ndarray of SeqSSRecord
        A 1D array of the precursor miRNA records sorted by indexes.
    """
    global _pre_mir_indexed_records
    if _pre_mir_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=PRECURSOR)
        records = {
            rec.id: rec for rec in iter_seq_ss(MIRNA_SEQ_SS_DIR, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _pre_mir_indexed_records = map_index(indexes)
    return _pre_mir_indexed_records


def index_mature_mirna(
    *, index_file: StrPath = INDEX_NPY
) -> NDArray[np.object_]:
    """Index the mature miRNA records.

    Parameters
    ----------
    index_file : str or os.PathLike, optional
        Path of the file containing the indexes. Default:
        `paths.INDEX_PY`.

    Returns
    -------
    numpy.ndarray of SeqSSRecord
        A 1D array of the mature miRNA records sorted by indexes.
    """
    global _mat_mir_indexed_records
    if _mat_mir_indexed_records is None:
        indexes = load_indexes(index_file=index_file, rna=MATURE)
        records = {
            rec.id: rec for rec in iter_fasta(MAT_MIRNA_FASTA, keep=indexes)
        }
        map_index = np.frompyfunc(records.__getitem__, 1, 1)
        _mat_mir_indexed_records = map_index(indexes)
    return _mat_mir_indexed_records


def interactions_from_index(
    *,
    interactions_csv: StrPath = PAIRS_FILE,
    index_file: StrPath = INDEX_NPY,
) -> NDArray[InteractionPair]:
    """Index the interaction pairs.

    Parameters
    ----------
    index_file : str or os.PathLike, optional
        Path of the file containing the indexes. Default:
        `paths.INDEX_PY`.

    Returns
    -------
    numpy.ndarray of InteractionPair
        A 1D array of the interaction pairs sorted by indexes.
    """
    global _indexed_interactions
    if _indexed_interactions is None:
        indexes = load_indexes(rna=INTERACTION)
        interactions = load_interactions(interactions_csv)
        it = (tuple(interaction.split("_")) for interaction in indexes)
        tuples = np.fromiter(
            it, dtype=dtype_interaction_pair, count=len(indexes)
        )
        mask = np.isin(interactions[PAIR_FIELDS], tuples)
        _indexed_interactions = interactions[mask]
    return _indexed_interactions


def dict_to_numpy(
    np_file: StrPath,
    feat_dict: dict[str, NDArray],
    *,
    index_file: StrPath = INDEX_NPY,
    rna: str = LNCRNA,
) -> None:
    """Convert a dictionary of numpy arrays to a indexed numpy array.

    Parameters
    ----------
    np_file : str or os.PathLike
        Path of the numpy file to save the indexed array.
    feat_dict : dict of str to numpy.ndarray
        Dictionary with the ID-array pairs to index.
    index_file : str or os.PathLike, optional
        Path of the index file to load. Default: `paths.INDEX_NPY`.
    rna : str, optional
        RNA type to use for indexing. Default: `LNCRNA`.
    """
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
    np.save(np_file, array, allow_pickle=False)
