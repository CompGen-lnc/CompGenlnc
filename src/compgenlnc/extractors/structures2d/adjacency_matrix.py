import math
import os
from collections.abc import Iterable, Iterator
from itertools import tee
from pathlib import Path

import numpy as np
import cv2

from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import iter_seq_ss, load_fasta
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SSLike


def get_adjacency_matrix(ss: SSLike) -> np.typing.NDArray[np.bool]:
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
    ss_list: Iterable[SSLike],
) -> Iterator[np.typing.NDArray[np.bool]]:
    for ss in ss_list:
        yield get_adjacency_matrix(ss)


def normalize_matrix_list(
    matrix_list: Iterable[np.typing.NDArray[np.bool]],
    normal_size: int | None = None,
) -> Iterator[np.typing.NDArray[np.float32]]:
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
    record_list: Iterable[SeqSSRecord],
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
) -> None:
    dict_file = Path(dict_file).resolve()
    dict_file.parent.mkdir(parents=True, exist_ok=True)

    id_vector = []
    ss_vector = []

    for record in record_list:
        id_vector.append(record.id)
        ss_vector.append(record.ss)

    mat_it = get_adjacency_matrix_list(ss_vector)
    normal_it = normalize_matrix_list(mat_it, normal_size)
    save_dict(dict_file, keys=id_vector, values=normal_it)


def gen_normalized_matrix_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
) -> None:
    filename = Path(filename).resolve()
    record_it = load_fasta(filename, mode="ss", keep=keep)
    gen_normalized_matrix_dict(record_it, dict_file, normal_size)


def gen_normalized_matrix_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
) -> None:
    folder = Path(folder).resolve()
    record_it = iter_seq_ss(folder, keep=keep)
    gen_normalized_matrix_dict(record_it, dict_file, normal_size)
