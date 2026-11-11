from __future__ import annotations

from collections.abc import Collection, Sequence
from typing import TYPE_CHECKING

import pandas as pd

from compgenlnc.fileman import load_fasta
from compgenlnc.structs import SeqSSRecord

if TYPE_CHECKING:
    from _typeshed import StrPath


def extract_precursor(
    mature: SeqSSRecord, /, precursors: Collection[SeqSSRecord]
) -> str:
    """Deduce the precursor miRNA for a mature miRNA.

    Parameters
    ----------
    mature : SeqSSRecord
       The mature miRNA record.
    precursors : collection of SeqSSRecord
        Set of precursor miRNA records.

    Returns
    -------
    str
        The precursor miRNA ID.
    """
    mat_id, mat_seq = mature.id, str(mature.seq)

    if not mat_id.endswith("p"):
        return mat_id
    base_precursor = mat_id.rsplit("-", 1)[0].lower()

    def precursor_matched(record):
        id_, seq = record.id.lower(), record.seq
        return base_precursor in id_ and mat_seq in seq

    match = next(filter(precursor_matched, precursors), None)
    return match.id if match else match


def gen_precursors_csv(
    csv: StrPath,
    matures: Sequence[SeqSSRecord],
    precursors: Collection[SeqSSRecord],
) -> None:
    """Generate a CSV file with miRNA associations.

    Parameters
    ----------
    csv : str or os.PathLike
        Path of the CSV file where associations will be saved.
    matures : sequence of SeqSSRecord
        Sequence of mature miRNA records.
    precursors : collection of SeqSSRecord
        Set of precursor miRNA records.
    """
    precursors = [extract_precursor(mir, precursors) for mir in matures]
    csv_dict = {
        "Mature": [rec.id for rec in matures],
        "Precursor": precursors,
    }
    df = pd.DataFrame.from_dict(csv_dict)
    df[df["Precursor"].notna()].to_csv(csv, index=False)


def gen_precursors_csv_from_files(
    csv: StrPath,
    mature_fasta: StrPath,
    precursor_fasta: StrPath,
) -> None:
    """Generate a CSV file with miRNA associations from FASTA files.

    Parameters
    ----------
    csv : str or os.PathLike
        Path of the CSV file where associations will be saved.
    mature_fasta : str or os.PathLike
        Path of the FASTA file containing mature miRNA sequences.
    precursor_fasta : str or os.PathLike
        Path of the FASTA file containing precursor miRNA sequences.
    """
    matures = load_fasta(mature_fasta)
    precursors = load_fasta(precursor_fasta)
    gen_precursors_csv(csv, matures, precursors)
