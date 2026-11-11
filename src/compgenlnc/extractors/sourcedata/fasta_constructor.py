from __future__ import annotations

import multiprocessing as mp
from collections.abc import Collection
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.params import SS_MODE
from compgenlnc.extractors.sourcedata.structure_predictor import (
    extract_2d_structure_identified,
)
from compgenlnc.fileman.fasta_manager import (
    iter_fasta,
    iter_seq_fasta,
    load_seq_ss,
    save_fasta,
)
from compgenlnc.fileman.log_manager import send_msg, track_progress

if TYPE_CHECKING:
    from _typeshed import StrPath


def filter_fasta(
    src: StrPath,
    dst: StrPath,
    /,
    filter: Collection[str],
) -> None:
    """Filter the RNA sequences from a FASTA file by the ID.

    Parameters
    ----------
    src : str or os.PathLike
        Path of the FASTA file.
    dst : str or os.PathLike
        Path of the filtered FASTA file.
    filter : collection of str
        A set with all substrings contained in the IDs of the desired
        sequences. The matched IDs must contain at least one substring
        from `filter`.
    """
    src = Path(src).resolve()
    dst = Path(dst).resolve()

    fasta = [record for record in iter_seq_fasta(src, filter)]
    fasta.sort(key=lambda record: record.id)

    dst.parent.mkdir(parents=True, exist_ok=True)
    save_fasta(dst, fasta, mode="seq")


def gen_fasta_2d(
    src: StrPath,
    dst: StrPath,
    /,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Generate a FASTA file containing the 2D structures calculated.

    Parameters
    ----------
    src : str or os.PathLike
        Path of the FASTA file with the RNA sequences which 2D structure
        will be calculated.
    dst : str or os.PathLike
        Path of the FASTA file where the 2D structures will be stored.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    send_msg("Generating 2D structures for sequences\n")

    if keep is None or not len(keep):
        total = sum(1 for _ in iter_fasta(src, mode="seq"))
    else:
        total = len(keep)
    dst = Path(dst).resolve()
    exclude = []

    if dst.exists():
        send_msg("  Analyzing existing FASTA...")
        exclude = [rec.id for rec in iter_fasta(dst, mode=SS_MODE)]
        send_msg("\r  Existing FASTA analyzed    \n")
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)

    num_processes = max(1, mp.cpu_count() - 2)
    with mp.Pool(processes=num_processes) as pool:
        ss_list = pool.imap_unordered(
            extract_2d_structure_identified,
            iter_fasta(src, mode="seq", keep=keep, exclude=exclude),
            chunksize=30,
        )
        initial = len(exclude)
        desc = "  Calculating 2D structures"
        ss_list = track_progress(ss_list, desc, total, "seq", initial)
        save_fasta(dst, ss_list, mode="ss")

    send_msg("  2D structure FASTA generated\n")


def seq_ss_to_fasta(
    src_dir: StrPath,
    dst_seq: StrPath,
    dst_ss: StrPath,
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Separate the sequence and 2D structure from seq+ss files.

    Parameters
    ----------
    src_dir : str or os.PathLike
        Path of the directory containing the seq+ss files.
    dst_seq : str or os.PathLike
        Path of the FASTA file where the sequences will be stored.
    dst_ss : str or os.PathLike
        Path of the FASTA file where the 2D structures will be stored.
    keep : collection of str, optional
        Set of IDs of the RNAs to save. Do not apply if `None` or empty.
        Default: `None`.
    """
    src_dir = Path(src_dir).resolve()
    dst_seq = Path(dst_seq).resolve()
    dst_ss = Path(dst_ss).resolve()

    dst_seq.parent.mkdir(parents=True, exist_ok=True)
    dst_ss.parent.mkdir(parents=True, exist_ok=True)
    records = load_seq_ss(src_dir, keep=keep)
    save_fasta(dst_seq, records, mode="seq")
    save_fasta(dst_ss, records, mode="ss")
