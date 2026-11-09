from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.params import LNCRNA, PRECURSOR
from compgenlnc.consts.paths import (
    DATASET_DIR,
    LNCRNA_DICTS_DIR,
    LNCRNA_LOOPS_DIR,
    LOOPS_DICT,
    MIRNA_LOOPS_DIR,
    PRE_MIRNA_DICTS_DIR,
)
from compgenlnc.dict_indexer import dict_to_numpy
from compgenlnc.extractors.structures2d.loops import gen_loops_dict
from compgenlnc.fileman.dict_manager import load_dict
from compgenlnc.structs.seq_ss_record import SeqSSRecord

if TYPE_CHECKING:
    from _typeshed import StrPath


def gen_ss_features(
    records: Sequence[SeqSSRecord],
    dicts_dir: StrPath,
    loops_dir: StrPath,
    *,
    rna: str,
) -> None:
    dicts_dir = Path(dicts_dir).resolve()
    dicts_dir.parent.mkdir(parents=True, exist_ok=True)

    loops_path = dicts_dir / LOOPS_DICT
    gen_loops_dict(loops_path, loops_dir, records)

    def path_to_npy(filename):
        return (DATASET_DIR / f"{rna}_{filename}").with_suffix(".npy")

    loops_npy = path_to_npy(LOOPS_DICT)
    dict_to_numpy(loops_npy, load_dict(loops_path), rna=rna)


def gen_lncrna_ss_features(records: Sequence[SeqSSRecord]) -> None:
    gen_ss_features(
        records,
        LNCRNA_DICTS_DIR,
        LNCRNA_LOOPS_DIR,
        rna=LNCRNA,
    )


def gen_mirna_ss_features(records: Sequence[SeqSSRecord]) -> None:
    gen_ss_features(
        records,
        PRE_MIRNA_DICTS_DIR,
        MIRNA_LOOPS_DIR,
        rna=PRECURSOR,
    )


def gen_ss_database(
    lnc_records: Sequence[SeqSSRecord],
    mir_records: Sequence[SeqSSRecord],
) -> None:
    gen_lncrna_ss_features(lnc_records)
    gen_mirna_ss_features(mir_records)
