from __future__ import annotations

from collections.abc import Collection
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.paths import (
    SEQ_SS_PREFIX,
)
from compgenlnc.fileman.fasta_manager import iter_fasta

if TYPE_CHECKING:
    from _typeshed import StrPath


def join_seq_ss(dir: StrPath, rna_id: str, seq: str, ss: str) -> None:
    """Save the sequence and 2D structure in the same file.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory where the seq+ss file will be stored.
    rna_id : str
        ID of the RNA.
    seq : str
        The RNA sequence on ACGT alphabet.
    ss : str
        The RNA 2d structure on dot-bracket notation.
    """
    dir = Path(dir).resolve()

    path = dir / f"{SEQ_SS_PREFIX}{rna_id}.dat"
    with open(path, "w") as out_file:
        out_file.write(seq + "\n")
        out_file.write(ss)


def fasta_to_seq_ss(
    src_seq: StrPath,
    src_ss: StrPath,
    dst_dir: StrPath,
    *,
    keep: Collection[str] | None = None,
) -> int:
    """Convert a FASTA file into seq+ss files for each RNA.

    Parameters
    ----------
    src_seq : str or os.PathLike
        Path of the FASTA file containing the RNA sequences.
    src_ss : str or os.PathLike
        Path of the FASTA file containing the RNA structures.
    dst_dir : str or os.PathLike
        Path of the directory where the seq+ss files will be saved.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    src_seq = Path(src_seq).resolve()
    src_ss = Path(src_ss).resolve()
    dst_dir = Path(dst_dir).resolve()

    seq_dict = {
        record.id: str(record.seq).upper().replace("U", "T")
        for record in iter_fasta(src_seq, mode="seq")
    }
    ss_dict = {
        record.id: str(record.seq) for record in iter_fasta(src_ss, mode="seq")
    }
    id_list = seq_dict.keys() & ss_dict.keys()
    if keep is not None:
        id_list &= set(keep)
    dst_dir.mkdir(parents=True, exist_ok=True)
    total = 0

    for id_ in id_list:
        try:
            seq = seq_dict[id_]
            ss = ss_dict[id_]
            join_seq_ss(dst_dir, id_, seq, ss)
            total += 1
        except Exception:
            pass

    return total
