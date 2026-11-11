from __future__ import annotations

import math
from collections.abc import Collection, Iterable, Iterator
from itertools import tee
from pathlib import Path
from typing import TYPE_CHECKING

import cv2
import numpy as np
from numpy.typing import NDArray

from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import iter_fasta, iter_seq_ss
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SSLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def get_adjacency_matrix(ss: SSLike, /) -> NDArray[np.bool]:
    """Get the adjacency matrix of a RNA structure.

    Parameters
    ----------
    ss : SSLike
        The 2D structure string on dot-bracket notation.

    Returns
    -------
    numpy.ndarray of bool
        A symmetric 2D matrix that indicates the paired nucleotides in
        the structure.
    """
    if isinstance(ss, SeqSSRecord):
        ss = ss.ss

    matrix = np.zeros((len(ss), len(ss)), np.bool)
    stack = []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            matrix[i, j] = matrix[j, i] = True

    return matrix


def get_adjacency_matrix_list(
    structures: Iterable[SSLike], /
) -> Iterator[NDArray[np.bool]]:
    """Get the adjacency matrices of a list of RNA structures.

    Parameters
    ----------
    structures : iterable object of SSLike
        Group of 2D structure strings on dot-bracket notation.

    Yields
    ------
    numpy.ndarray of bool
        A symmetric 2D matrix that indicates the paired nucleotides in
        the structure.
    """
    for ss in structures:
        yield get_adjacency_matrix(ss)


def normalize_matrix_list(
    matrices: Iterable[NDArray[np.bool]],
    /,
    normal_size: int | None = None,
) -> Iterator[NDArray[np.float32]]:
    """Normalize a list of adjacency matrix.

    Parameters
    ----------
    matrices : iterable object of numpy.ndarray of bool
        Group of adjacency matrices.
    normal_size : int, optional
        Size used to resize the adjacency matrices. If `None`, it will
        use the mean length of the matrices. Default: `None`.

    Yields
    ------
    numpy.ndarray of bool
        A normalized 2D matrix rescaled with bilinear interpolation.
    """
    if normal_size is None:
        matrices, matrix_copy = tee(matrices)
        normal_size = math.ceil(np.mean([len(mat) for mat in matrix_copy]))

    for mat in matrices:
        yield cv2.resize(
            mat.astype(np.uint8),
            (normal_size, normal_size),
            interpolation=cv2.INTER_LINEAR,
        )


def gen_normalized_matrix_dict(
    dict_file: StrPath,
    records: Iterable[SeqSSRecord],
    normal_size: int | None = None,
) -> None:
    """Generate a dictionary of normalized matrices.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    records : iterable object of SeqSSRecord
        Set containing all the 2D structure records to save into the
        dictionary.
    normal_size : int, optional
        Size used to resize the adjacency matrices. If `None`, it will
        use the mean length of the matrices. Default: `None`.
    """
    dict_file = Path(dict_file).resolve()
    dict_file.parent.mkdir(parents=True, exist_ok=True)

    id_vector = []
    ss_vector = []

    for record in records:
        id_vector.append(record.id)
        ss_vector.append(record.ss)

    if normal_size is None:
        normal_size = math.ceil(np.mean([len(ss) for ss in ss_vector]))
    mat_it = get_adjacency_matrix_list(ss_vector)
    normal_it = normalize_matrix_list(mat_it, normal_size)
    save_dict(dict_file, keys=id_vector, values=normal_it, overwrite=True)


def gen_normalized_matrix_dict_from_fasta(
    dict_file: StrPath,
    fasta_file: StrPath,
    normal_size: int | None = None,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of normalized matrices from a FASTA file.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    fasta_file : str or os.PathLike
        Path of the FASTA file with the 2D structures.
    normal_size : int, optional
        Size used to resize the adjacency matrices. If `None`, it will
        use the mean length of the matrices. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    records = iter_fasta(fasta_file, mode="ss", keep=keep)
    gen_normalized_matrix_dict(dict_file, records, normal_size)


def gen_normalized_matrix_dict_from_dir(
    dict_file: StrPath,
    seq_ss_dir: StrPath,
    normal_size: int | None = None,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of normalized matrices from seq+ss files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    normal_size : int, optional
        Size used to resize the adjacency matrices. If `None`, it will
        use the mean length of the matrices. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    records = iter_seq_ss(seq_ss_dir, keep=keep)
    gen_normalized_matrix_dict(dict_file, records, normal_size)
