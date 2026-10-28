import os
from collections.abc import Sequence
from pathlib import Path

from compgenlnc.consts.params import LNCRNA, PRECURSOR
from compgenlnc.consts.paths import (
    DATASET_FOLDER,
    LNCRNA_DICTS_FOLDER,
    LNCRNA_LOOPS_FOLDER,
    LOOPS_DICT,
    MIRNA_LOOPS_FOLDER,
    PRE_MIRNA_DICTS_FOLDER,
)
from compgenlnc.dict_indexer import dict_to_numpy
from compgenlnc.extractors.structures2d.loops import gen_loops_dict
from compgenlnc.fileman.dict_manager import load_dict
from compgenlnc.structs.seq_ss_record import SeqSSRecord


def gen_ss_features(
    record_list: Sequence[SeqSSRecord],
    dicts_folder: str | os.PathLike,
    loops_folder: str | os.PathLike,
    *,
    rna: str,
) -> None:
    dicts_folder = Path(dicts_folder).resolve()
    dicts_folder.parent.mkdir(parents=True, exist_ok=True)

    loops_path = dicts_folder / LOOPS_DICT
    gen_loops_dict(loops_path, record_list, loops_folder)

    def path_to_npy(filename):
        return (DATASET_FOLDER / f"{rna}_{filename}").with_suffix(".npy")

    loops_npy = path_to_npy(LOOPS_DICT)
    dict_to_numpy(loops_npy, load_dict(loops_path), rna=rna)


def gen_lncrna_ss_features(record_list: Sequence[SeqSSRecord]) -> None:
    gen_ss_features(
        record_list,
        LNCRNA_DICTS_FOLDER,
        LNCRNA_LOOPS_FOLDER,
        rna=LNCRNA,
    )


def gen_mirna_ss_features(record_list: Sequence[SeqSSRecord]) -> None:
    gen_ss_features(
        record_list,
        PRE_MIRNA_DICTS_FOLDER,
        MIRNA_LOOPS_FOLDER,
        rna=PRECURSOR,
    )


def gen_ss_database(
    lnc_records: Sequence[SeqSSRecord],
    mir_records: Sequence[SeqSSRecord],
) -> None:
    gen_lncrna_ss_features(lnc_records)
    gen_mirna_ss_features(mir_records)
