from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from compgenlnc.consts.fields import LNCRNA_FIELD, MIRNA_FIELD
from compgenlnc.consts.params import INTERACTION, LNCRNA, MIRNA
from compgenlnc.consts.paths import (
    BINDING_PREFIX,
    CTD_DICT,
    DATASET_DIR,
    DOC2VEC_DICT,
    KMER_DICT,
    LNCRNA_DICTS_DIR,
    LNCRNA_INTERACTION_DOV2VEC_MODEL,
    LNCRNA_LOOPS_DIR,
    LOOPS_DICT,
    MIRNA_DICTS_DIR,
    MIRNA_INTERACTION_DOV2VEC_MODEL,
    MIRNA_LOOPS_DIR,
    NORMAL_SEQ_DICT,
    ROLE2VEC_DICT,
)
from compgenlnc.dict_indexer import dict_to_numpy
from compgenlnc.extractors.interactions import gen_substructure_loops_dict
from compgenlnc.extractors.sequences.ctd import gen_ctd_dict
from compgenlnc.extractors.sequences.doc2vec import gen_doc2vec_dict
from compgenlnc.extractors.sequences.kmer import gen_kmer_dict
from compgenlnc.extractors.sequences.normal_sequence import (
    gen_normalized_sequence_dict,
)
from compgenlnc.extractors.sequences.role2vec import (
    gen_role2vec_dict_from_files,
)
from compgenlnc.fileman.dict_manager import load_dict
from compgenlnc.fileman.loops_manager import read_2d_structure_loops
from compgenlnc.fileman.miranda_manager import get_subrecord_list
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.structs.loop_counter import LoopCounter
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import InteractionPair

if TYPE_CHECKING:
    from _typeshed import StrPath


def gen_interaction_features(
    records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    dicts_dir: StrPath,
    loops_counters: Sequence[LoopCounter],
    binding_zones: Sequence[BindingZone],
    k: int,
    dv2_model_file: StrPath,
    *,
    rna: str,
) -> None:
    dicts_dir = Path(dicts_dir).resolve()
    dicts_dir.parent.mkdir(parents=True, exist_ok=True)

    kmer_path = dicts_dir / f"{BINDING_PREFIX}{KMER_DICT}"
    ctd_path = dicts_dir / f"{BINDING_PREFIX}{CTD_DICT}"
    d2v_path = dicts_dir / f"{BINDING_PREFIX}{DOC2VEC_DICT}"
    r2v_path = dicts_dir / f"{BINDING_PREFIX}{ROLE2VEC_DICT}"
    normal_path = dicts_dir / f"{BINDING_PREFIX}{NORMAL_SEQ_DICT}"
    loops_path = dicts_dir / f"{BINDING_PREFIX}{LOOPS_DICT}"

    gen_kmer_dict(kmer_path, records, k)
    gen_ctd_dict(ctd_path, records)
    gen_doc2vec_dict(d2v_path, records, model_file=dv2_model_file)
    gen_role2vec_dict_from_files(r2v_path, kmer_path, ctd_path, d2v_path)
    gen_normalized_sequence_dict(normal_path, records)
    gen_substructure_loops_dict(
        loops_path, interactions, loops_counters, binding_zones, molecule=rna
    )

    def path_to_npy(filename):
        path = DATASET_DIR / f"binding_{rna}_{filename}"
        return path.with_suffix(".npy")

    kmer_npy = path_to_npy(KMER_DICT)
    ctd_npy = path_to_npy(CTD_DICT)
    d2v_npy = path_to_npy(DOC2VEC_DICT)
    r2v_npy = path_to_npy(ROLE2VEC_DICT)
    normal_npy = path_to_npy(NORMAL_SEQ_DICT)
    loops_npy = path_to_npy(LOOPS_DICT)
    dict_to_numpy(kmer_npy, load_dict(kmer_path), rna=INTERACTION)
    dict_to_numpy(ctd_npy, load_dict(ctd_path), rna=INTERACTION)
    dict_to_numpy(d2v_npy, load_dict(d2v_path), rna=INTERACTION)
    dict_to_numpy(r2v_npy, load_dict(r2v_path), rna=INTERACTION)
    dict_to_numpy(normal_npy, load_dict(normal_path), rna=INTERACTION)
    dict_to_numpy(loops_npy, load_dict(loops_path), rna=INTERACTION)


def gen_lncrna_interaction_features(
    records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    loops_counters: Sequence[LoopCounter],
    binding_zones: Sequence[BindingZone],
) -> None:
    gen_interaction_features(
        records,
        interactions,
        LNCRNA_DICTS_DIR,
        loops_counters,
        binding_zones,
        4,
        LNCRNA_INTERACTION_DOV2VEC_MODEL,
        rna=LNCRNA,
    )


def gen_mirna_interaction_features(
    records: Sequence[SeqSSRecord],
    pre_mir_records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    loops_counters: Sequence[LoopCounter],
    binding_zones: Sequence[BindingZone],
) -> None:
    def adjust_binding_list():
        for record, precursor, binding_zone in zip(
            records, pre_mir_records, binding_zones
        ):
            # print(precursor.id, record.id)
            index = str(precursor.seq).index(str(record.seq))
            yield binding_zone.adjust_mirna_range(index)

    for binding in adjust_binding_list():
        pass

    gen_interaction_features(
        records,
        interactions,
        MIRNA_DICTS_DIR,
        loops_counters,
        (binding for binding in adjust_binding_list()),
        # binding_zones,
        3,
        MIRNA_INTERACTION_DOV2VEC_MODEL,
        rna=MIRNA,
    )


def gen_interaction_database(
    lnc_records: Sequence[SeqSSRecord],
    mir_records: Sequence[SeqSSRecord],
    pre_mir_records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    precursor_map: dict[str, str],
    binding_zones: Sequence[BindingZone],
) -> None:
    def subrecords(records, rna) -> Sequence[SeqSSRecord]:
        count = len(records)
        subrecord_list = np.fromiter(
            get_subrecord_list(records, interactions, binding_zones, rna=rna),
            dtype=object,
            count=count,
        )
        return subrecord_list

    gen_lncrna_interaction_features(
        subrecords(lnc_records, LNCRNA),
        interactions,
        (
            read_2d_structure_loops(LNCRNA_LOOPS_DIR, lnc)
            for lnc in interactions[LNCRNA_FIELD]
        ),
        binding_zones,
    )
    gen_mirna_interaction_features(
        subrecords(mir_records, MIRNA),
        pre_mir_records,
        interactions,
        (
            read_2d_structure_loops(MIRNA_LOOPS_DIR, precursor_map[mir])
            for mir in interactions[MIRNA_FIELD]
        ),
        binding_zones,
    )
