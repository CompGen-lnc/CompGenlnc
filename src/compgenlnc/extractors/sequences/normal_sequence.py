import math
import os
from collections.abc import Iterable, Iterator
from itertools import tee
from pathlib import Path

import numpy as np

from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.fasta_manager import iter_seq_ss, load_fasta


def recode_sequence(seq: SeqLike) -> np.typing.NDArray[np.uint8]:
    if isinstance(seq, SeqSSRecord):
        seq = seq.seq

    recode_map = {"A": 1, "C": 2, "G": 3, "T": 4, "U": 4}
    recoded = np.fromiter(
        (recode_map[char] if char in recode_map.keys() else 5 for char in seq),
        np.uint8,
        len(seq),
    )

    return recoded


def recode_sequence_list(
    seq_list: Iterable[SeqLike],
) -> Iterator[np.typing.NDArray[np.uint8]]:
    for seq in seq_list:
        yield recode_sequence(seq)


def normalize_sequence_list(
    seq_list: Iterable[np.typing.NDArray[np.uint8]],
    normal_size: int | None = None,
) -> np.typing.NDArray[np.float32]:
    if normal_size is None:
        seq_list, seq_copy = tee(seq_list)
        normal_size = math.ceil(np.mean([len(seq) for seq in seq_copy]))

    matrix = np.array(
        [
            np.pad(
                sequence[:normal_size],
                (0, max(0, normal_size - len(sequence))),
            )
            for sequence in seq_list
        ],
        np.float32,
    )

    average = matrix.mean(0)
    deviation = matrix.std(0)
    deviation[deviation < 1e-10] = 1

    return (matrix - average) / deviation


def gen_normalized_sequence_dict(
    record_list: Iterable[SeqSSRecord],
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
) -> None:
    dict_file = Path(dict_file).resolve()
    dict_file.parent.mkdir(parents=True, exist_ok=True)

    id_vector = []
    seq_vector = []

    for record in record_list:
        id_vector.append(record.id)
        seq_vector.append(recode_sequence(record.seq))

    seq_dict = {
        pair[0]: pair[1]
        for pair in zip(
            id_vector, normalize_sequence_list(seq_vector, normal_size)
        )
    }
    save_dict(dict_file, seq_dict)


def gen_normalized_sequence_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    filename = Path(filename).resolve()
    dict_file = Path(dict_file).resolve()

    dict_file.parent.mkdir(parents=True, exist_ok=True)

    id_vector = []
    seq_vector = []

    for record in load_fasta(filename, mode="seq", keep=keep):
        id_vector.append(record.id)
        seq_vector.append(recode_sequence(record.seq))

    seq_dict = {
        pair[0]: pair[1]
        for pair in zip(
            id_vector, normalize_sequence_list(seq_vector, normal_size)
        )
    }
    save_dict(dict_file, seq_dict)


def gen_normalized_sequence_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_file.parent.mkdir(parents=True, exist_ok=True)

    id_vector = []
    seq_vector = []

    for record in iter_seq_ss(folder, keep):
        id_vector.append(record.id)
        seq_vector.append(recode_sequence(record.seq))

    seq_dict = {
        pair[0]: pair[1]
        for pair in zip(
            id_vector, normalize_sequence_list(seq_vector, normal_size)
        )
    }
    save_dict(dict_file, seq_dict)
