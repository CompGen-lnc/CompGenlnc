from collections import Counter
from itertools import product
from pathlib import Path
from typing import Iterable
import os

from Bio.SeqRecord import SeqRecord
import numpy as np

from compgenlnc.config.paths import SEQ_SS_PREFIX
from compgenlnc.utils.dict_manager import save_dict
from compgenlnc.utils.fasta_manager import read_sequence, iter_seq_ss


# Array with all possible combinations for RNA k-mers with k from 1 to 4
n_combinations = [
    [''.join(p) for p in product('ATCG', repeat=n)]
    for n in range(1, 5)
]


def get_kmer(seq: str | SeqRecord, k: int) -> np.typing.NDArray[np.float64]:
    """Get the k-mers of a RNA sequences from 1 to k from the given
    sequence.

    Args:
        seq: The RNA sequence string on ATCG alphabet.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.

    Returns:
        NDArray: A concatenated numpy array with the k-mers of the
            sequence.
    """
    # Get the kmers of length n for the sequence
    if isinstance(seq, SeqRecord):
        seq = str(seq.seq)

    seq = seq.upper().replace('U', 'T')
    def get_nmer(seq: str, n: int) -> np.typing.NDArray[np.float64]:
        if not seq or len(seq) < n:
            return np.array([0] * 4 ** n)
        kmers = Counter(
            (seq[i : i + n] for i in range(len(seq) - n + 1))
        )
        combinations = n_combinations[n - 1]
        return np.array([
            kmers[sub_s] / kmers.total()
            for sub_s in combinations
        ])
    return np.concatenate([get_nmer(seq, n) for n in range(1, k + 1)])


def get_kmer_by_name(
        id_: str,
        folder: str | os.PathLike,
        k: int,
        /,
        mature_only: bool = False
) -> np.typing.NDArray[np.float64]:
    """Get the k-mers of a RNA sequences from 1 to k from the seq+ss
    file of the given molecule.

    Args:
        id_: Name of the RNA sequence.
        folder: Path of the folder where the seq+ss file is.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        mature_only: Only saves the mature section of the sequence if
            True.

    Returns:
        NDArray: A concatenated numpy array with the k-mers of the
            sequence.
    """
    folder = Path(folder).resolve()
    
    filename = folder / f'{SEQ_SS_PREFIX}{id_}.dat'
    seq = read_sequence(filename, mature_only=mature_only)
    return get_kmer(seq, k)


def gen_kmer_dict(
        record_list: Iterable[SeqRecord],
        dict_file: str | os.PathLike,
        k: int,
        /,
        mature_only: bool = False
) -> None:
    """Get the k-mers of a RNA sequences from 1 to k from the seq+ss
    file of the given molecule.

    Args:
        record_list: Iterable that contains all the sequence records to
            save into the dictionary.
        dict_file: Path of the file where the dictionary will be saved.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    dict_file = Path(dict_file).resolve()
    
    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    kmer_dict = {record.id: get_kmer(record.seq, k) for record in record_list}
    save_dict(dict_file, kmer_dict)


def gen_kmer_dict_from_folder(
        folder: str | os.PathLike,
        dict_file: str | os.PathLike,
        k: int,
        /,
        keep: Iterable[str] | None = None,
        mature_only: bool = False
) -> None:
    """Get the k-mers of a RNA sequences from 1 to k from the seq+ss
    file of the given molecule.

    Args:
        folder: Path of the folder where the seq+ss file is.
        dict_file: Path of the file where the dictionary will be saved.
        k: Maximum length for k-mers. It calculates k-mers with length
            from 1 to k, included.
        keep: List of the RNA names to save. Save all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    folder = Path(folder).resolve()
    dict_file = Path(dict_file).resolve()

    dict_folder = dict_file.parent
    dict_folder.mkdir(parents=True, exist_ok=True)

    kmer_dict = {
        record.id: get_kmer(record.seq, k)
        for record in iter_seq_ss(folder, keep=keep, mature_only=mature_only)
    }
    save_dict(dict_file, kmer_dict)
