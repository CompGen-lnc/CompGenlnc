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
    matures = load_fasta(mature_fasta)
    precursors = load_fasta(precursor_fasta)
    gen_precursors_csv(csv, matures, precursors)
