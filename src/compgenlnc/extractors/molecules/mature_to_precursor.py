import os
from collections.abc import Iterable, Sequence

import pandas as pd

from compgenlnc.fileman import load_fasta
from compgenlnc.structs import SeqSSRecord


def extract_precursor(
    mir: str | SeqSSRecord, precursor_list: Iterable[str]
) -> str:
    if isinstance(mir, SeqSSRecord):
        mir = mir.id

    if not mir.endswith("p"):
        return mir

    last_hyphen = mir.rfind("-")
    base_precursor = mir[:last_hyphen]
    return (id_ for id_ in precursor_list if base_precursor in id_)


def gen_precursors_csv(
    csv: str | os.PathLike,
    mir_it: Iterable[str | SeqSSRecord],
    precursor_it: Iterable[str | SeqSSRecord],
) -> None:
    mir_list = [
        mir.id if isinstance(mir, SeqSSRecord) else mir for mir in mir_it
    ]
    precursor_list = [
        pre_.id if isinstance(pre_, SeqSSRecord) else pre_
        for pre_ in precursor_it
    ]
    precursors = (extract_precursor(mir, precursor_list) for mir in mir_list)
    csv_dict = {"Mature": mir_list, "Precursor": precursors}
    df = pd.DataFrame.from_dict(csv_dict)
    df.to_csv(csv)


def gen_precursors_csv_from_files(
    csv: str | os.PathLike,
    mature_fasta: str | os.PathLike,
    precursor_fasta: str | os.PathLike,
) -> None:
    mature_list = load_fasta(mature_fasta)
    precursor_list = load_fasta(precursor_fasta)
    gen_precursors_csv(csv, mature_list, precursor_list)
