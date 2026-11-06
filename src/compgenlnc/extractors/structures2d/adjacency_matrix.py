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
from compgenlnc.fileman.fasta_manager import iter_seq_ss, iter_fasta
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SSLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def get_adjacency_matrix(ss: SSLike, /) -> NDArray[np.bool]:
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
    ss_list: Iterable[SSLike], /
) -> Iterator[NDArray[np.bool]]:
    for ss in ss_list:
        yield get_adjacency_matrix(ss)


def normalize_matrix_list(
    matrix_list: Iterable[NDArray[np.bool]],
    /,
    normal_size: int | None = None,
) -> Iterator[NDArray[np.float32]]:
    if normal_size is None:
        matrix_list, matrix_copy = tee(matrix_list)
        normal_size = math.ceil(np.mean([len(mat) for mat in matrix_copy]))

    for mat in matrix_list:
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
    records = iter_fasta(fasta_file, mode="ss", keep=keep)
    gen_normalized_matrix_dict(dict_file, records, normal_size)


def gen_normalized_matrix_dict_from_folder(
    dict_file: StrPath,
    seq_ss_folder: StrPath,
    normal_size: int | None = None,
    *,
    keep: Collection[str] | None = None,
) -> None:
    records = iter_seq_ss(seq_ss_folder, keep=keep)
    gen_normalized_matrix_dict(dict_file, records, normal_size)
