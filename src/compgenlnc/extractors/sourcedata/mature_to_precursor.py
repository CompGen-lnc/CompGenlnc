import os
from collections.abc import Iterable, Sequence

import pandas as pd

from compgenlnc.fileman import load_fasta
from compgenlnc.structs import SeqSSRecord


def extract_precursor(
    mir: SeqSSRecord, precursor_list: Iterable[SeqSSRecord]
) -> str:
    mat_id, mat_seq = mir.id, str(mir.seq)

    if not mat_id.endswith("p"):
        return mat_id
    base_precursor = mat_id.rsplit("-", 1)[0].lower()

    def precursor_matched(record):
        id_, seq = record.id.lower(), record.seq
        return base_precursor in id_ and mat_seq in seq

    match = next(filter(precursor_matched, precursor_list), None)
    return match.id if match else match


def gen_precursors_csv(
    csv: str | os.PathLike,
    mir_it: Iterable[SeqSSRecord],
    precursor_it: Iterable[SeqSSRecord],
) -> None:
    mir_list = list(mir_it)
    precursor_list = list(precursor_it)
    precursors = [extract_precursor(mir, precursor_list) for mir in mir_list]
    csv_dict = {
        "Mature": [rec.id for rec in mir_list],
        "Precursor": precursors,
    }
    df = pd.DataFrame.from_dict(csv_dict)
    df[df["Precursor"].notna()].to_csv(csv, index=False)


def gen_precursors_csv_from_files(
    csv: str | os.PathLike,
    mature_fasta: str | os.PathLike,
    precursor_fasta: str | os.PathLike,
) -> None:
    mature_list = load_fasta(mature_fasta)
    precursor_list = load_fasta(precursor_fasta)
    gen_precursors_csv(csv, mature_list, precursor_list)
