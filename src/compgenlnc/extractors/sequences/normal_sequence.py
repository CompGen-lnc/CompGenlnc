from __future__ import annotations

import math
from collections.abc import Collection, Iterable, Iterator, Sequence
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from compgenlnc.consts.params import SEQ_MODE
from compgenlnc.fileman.dict_manager import save_dict
from compgenlnc.fileman.fasta_manager import iter_fasta, load_seq_ss
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.molecules import SeqLike

if TYPE_CHECKING:
    from _typeshed import StrPath


def recode_sequence(seq: SeqLike, /) -> NDArray[np.uint8]:
    """Recode a sequence into a numerical vector.

    Parameters
    ----------
    seq : SeqLike
        The RNA sequence string on AUCG alphabet.

    Returns
    -------
    numpy.ndarray of uint8
        A 1D array with the recoded sequence.
    """
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
    sequences: Iterable[SeqLike], /
) -> Iterator[NDArray[np.uint8]]:
    """Recode a list of sequences into numerical vectors.

    Parameters
    ----------
    sequences : iterable object of SeqLike
        Group of RNA sequence strings on AUCG alphabet.

    Yields
    ------
    numpy.ndarray of uint8
        A 1D array with the recoded sequence.
    """
    for seq in sequences:
        yield recode_sequence(seq)


def normalize_sequence_list(
    sequences: Sequence[NDArray[np.uint8]],
    /,
    normal_size: int | None = None,
) -> NDArray[np.float32]:
    """Normalize a list of recoded sequences.

    Parameters
    ----------
    sequences : sequence of numpy.ndarray of uint8
        List of recoded sequences.
    normal_size : int, optional
        Size used to resize the recoded sequences. If `None`, it will
        use the mean length of the sequences. Default: `None`.

    Returns
    -------
    numpy.ndarray of float32
        A 2D matrix with the normalized sequences as rows.
    """
    if normal_size is None:
        normal_size = math.ceil(np.mean([len(seq) for seq in sequences]))

    matrix = np.array(
        [
            np.pad(
                sequence[:normal_size],
                (0, max(0, normal_size - len(sequence))),
            )
            for sequence in sequences
        ],
        np.float32,
    )

    average = matrix.mean(0)
    deviation = matrix.std(0)
    deviation[deviation < 1e-10] = 1

    return (matrix - average) / deviation


def gen_normalized_sequence_dict(
    dict_file: StrPath,
    records: Iterable[SeqSSRecord],
    normal_size: int | None = None,
) -> None:
    """Generate a dictionary of normalized sequences.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    records : iterable object of SeqSSRecord
        Set containing all the sequence records to save into the
        dictionary.
    normal_size : int, optional
        Size used to resize the recoded sequences. If `None`, it will
        use the mean length of the sequences. Default: `None`.
    """
    send_msg("Generating normalized sequences dictionary\n")
    id_vector = []
    seq_vector = []

    desc = "  Recoding sequences"
    records = track_progress(records, desc=desc)
    for record in records:
        id_vector.append(record.id)
        seq_vector.append(recode_sequence(record.seq))

    desc = "  Normalizing sequences"
    pairs = zip(id_vector, normalize_sequence_list(seq_vector, normal_size))
    pairs = track_progress(pairs, desc, unit="seq")
    seq_dict = {pair[0]: pair[1] for pair in pairs}
    save_dict(dict_file, seq_dict, overwrite=True)
    send_msg("  Normalized sequences dictionary generated\n")


def gen_normalized_sequence_dict_from_fasta(
    dict_file: StrPath,
    fasta_file: StrPath,
    normal_size: int | None = None,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of normalized sequences from a FASTA file.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    fasta_file : str or os.PathLike
        Path of the FASTA file with the sequences.
    normal_size : int, optional
        Size used to resize the recoded sequences. If `None`, it will
        use the mean length of the sequences. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    records = iter_fasta(fasta_file, mode=SEQ_MODE, keep=keep)
    gen_normalized_sequence_dict(dict_file, records, normal_size)


def gen_normalized_sequence_dict_from_dir(
    dict_file: StrPath,
    seq_ss_dir: StrPath,
    normal_size: int | None = None,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Generate a dictionary of normalized sequences from seq+ss files.

    Parameters
    ----------
    dict_file : str or os.PathLike
        Path of the file where the dictionary will be saved.
    seq_ss_dir : str or os.PathLike
        Path of the directory where the seq+ss files are.
    normal_size : int, optional
        Size used to resize the recoded sequences. If `None`, it will
        use the mean length of the sequences. Default: `None`.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    records = load_seq_ss(seq_ss_dir, keep=keep)
    gen_normalized_sequence_dict(dict_file, records, normal_size)
