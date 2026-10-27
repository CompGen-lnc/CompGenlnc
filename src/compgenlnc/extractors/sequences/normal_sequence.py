import math
import os
from collections.abc import Iterable, Iterator
from itertools import tee

import numpy as np

from compgenlnc.consts.params import SEQ_MODE
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import iter_seq_ss, load_fasta
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike


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
    send_msg("Generating normalized sequences dictionary\n")
    id_vector = []
    seq_vector = []

    desc = "  Recoding sequences"
    record_list = track_progress(record_list, desc=desc)
    for record in record_list:
        id_vector.append(record.id)
        seq_vector.append(recode_sequence(record.seq))

    desc = "  Normalizing sequences"
    pairs = zip(id_vector, normalize_sequence_list(seq_vector, normal_size))
    pairs = track_progress(pairs, desc=desc)
    seq_dict = {pair[0]: pair[1] for pair in pairs}
    save_dict(dict_file, seq_dict, overwrite=True)
    send_msg("  Normalized sequences dictionary generated\n")


def gen_normalized_sequence_dict_from_fasta(
    filename: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    record_list = load_fasta(filename, mode=SEQ_MODE, keep=keep)
    gen_normalized_sequence_dict(record_list, dict_file, normal_size)


def gen_normalized_sequence_dict_from_folder(
    folder: str | os.PathLike,
    dict_file: str | os.PathLike,
    normal_size: int | None = None,
    /,
    keep: Iterable[str] | None = None,
    mature_only: bool = False,
) -> None:
    record_list = iter_seq_ss(folder, keep=keep)
    gen_normalized_sequence_dict(record_list, dict_file, normal_size)
