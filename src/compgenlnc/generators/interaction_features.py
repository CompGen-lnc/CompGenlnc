import os
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

from compgenlnc.consts.paths import (
    BINDING_PREFIX,
    CTD_DICT,
    DATASET_FOLDER,
    DOC2VEC_DICT,
    KMER_DICT,
    LNCRNA_DICTS_FOLDER,
    LNCRNA_INTERACTION_DOV2VEC_MODEL,
    LNCRNA_LOOPS_FOLDER,
    LOOPS_DICT,
    LOOPS_PREFIX,
    MIRNA_DICTS_FOLDER,
    MIRNA_INTERACTION_DOV2VEC_MODEL,
    MIRNA_LOOPS_FOLDER,
    NORMAL_SEQ_DICT,
    ROLE2VEC_DICT,
)
from compgenlnc.consts.fields import LNCRNA_FIELD, MIRNA_FIELD
from compgenlnc.consts.params import INTERACTION, LNCRNA, MIRNA
from compgenlnc.dict_indexer import dict_to_numpy
from compgenlnc.extractors.sequences.ctd import gen_ctd_dict
from compgenlnc.extractors.sequences.doc2vec import gen_doc2vec_dict
from compgenlnc.extractors.sequences.kmer import gen_kmer_dict
from compgenlnc.extractors.sequences.normal_sequence import (
    gen_normalized_sequence_dict,
)
from compgenlnc.extractors.sequences.role2vec import (
    gen_role2vec_dict_from_files,
)
from compgenlnc.extractors.interactions import gen_substructure_loops_dict
from compgenlnc.fileman.dict_manager import load_dict
from compgenlnc.fileman.loops_manager import read_2d_structure_loops
from compgenlnc.fileman.miranda_manager import get_subrecord_list
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.structs.loop_counter import LoopCounter
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import InteractionPair


def gen_interaction_features(
    record_list: Sequence[SeqSSRecord],
    interaction_list: Sequence[InteractionPair],
    dicts_folder: str | os.PathLike,
    loops_list: Iterable[LoopCounter],
    binding_list: Iterable[BindingZone],
    k: int,
    dv2_model_file: str | os.PathLike,
    *,
    rna: str,
) -> None:
    dicts_folder = Path(dicts_folder).resolve()
    dicts_folder.parent.mkdir(parents=True, exist_ok=True)

    kmer_path = dicts_folder / f"{BINDING_PREFIX}{KMER_DICT}"
    ctd_path = dicts_folder / f"{BINDING_PREFIX}{CTD_DICT}"
    d2v_path = dicts_folder / f"{BINDING_PREFIX}{DOC2VEC_DICT}"
    r2v_path = dicts_folder / f"{BINDING_PREFIX}{ROLE2VEC_DICT}"
    normal_path = dicts_folder / f"{BINDING_PREFIX}{NORMAL_SEQ_DICT}"
    loops_path = dicts_folder / f"{BINDING_PREFIX}{LOOPS_DICT}"

    gen_kmer_dict(record_list, kmer_path, k)
    gen_ctd_dict(record_list, ctd_path)
    gen_doc2vec_dict(record_list, d2v_path, model_file=dv2_model_file)
    gen_role2vec_dict_from_files(kmer_path, ctd_path, d2v_path, r2v_path)
    gen_normalized_sequence_dict(record_list, normal_path)
    gen_substructure_loops_dict(
        loops_path, interaction_list, loops_list, binding_list, molecule=rna
    )

    def path_to_npy(filename):
        path = DATASET_FOLDER / f"binding_{rna}_{filename}"
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
    record_list: Sequence[SeqSSRecord],
    interaction_list: Sequence[InteractionPair],
    loops_list: Iterable[LoopCounter],
    binding_list: Iterable[BindingZone],
) -> None:
    gen_interaction_features(
        record_list,
        interaction_list,
        LNCRNA_DICTS_FOLDER,
        loops_list,
        binding_list,
        4,
        LNCRNA_INTERACTION_DOV2VEC_MODEL,
        rna=LNCRNA,
    )


def gen_mirna_interaction_features(
    record_list: Sequence[SeqSSRecord],
    pre_mir_records: Sequence[SeqSSRecord],
    interaction_list: Sequence[InteractionPair],
    loops_list: Iterable[LoopCounter],
    binding_list: Iterable[BindingZone],
) -> None:
    def adjust_binding_list():
        for record, precursor, binding_zone in zip(
            record_list, pre_mir_records, binding_list
        ):
            # print(precursor.id, record.id)
            index = str(precursor.seq).index(str(record.seq))
            yield binding_zone.adjust_mirna_range(index)

    for binding in adjust_binding_list():
        pass

    gen_interaction_features(
        record_list,
        interaction_list,
        MIRNA_DICTS_FOLDER,
        loops_list,
        (binding for binding in adjust_binding_list()),
        # binding_list,
        3,
        MIRNA_INTERACTION_DOV2VEC_MODEL,
        rna=MIRNA,
    )


def gen_interaction_database(
    lnc_records: Sequence[SeqSSRecord],
    mir_records: Sequence[SeqSSRecord],
    pre_mir_records: Sequence[SeqSSRecord],
    interaction_list: Sequence[InteractionPair],
    precursor_map: dict[str, str],
    binding_list: Sequence[BindingZone],
) -> None:
    def subrecords(record_list, rna) -> Sequence[SeqSSRecord]:
        count = len(record_list)
        subrecord_list = np.fromiter(
            get_subrecord_list(
                record_list, interaction_list, binding_list, rna=rna
            ),
            dtype=object,
            count=count,
        )
        return subrecord_list

    gen_lncrna_interaction_features(
        subrecords(lnc_records, LNCRNA),
        interaction_list,
        (
            read_2d_structure_loops(LNCRNA_LOOPS_FOLDER, lnc)
            for lnc in interaction_list[LNCRNA_FIELD]
        ),
        binding_list,
    )
    gen_mirna_interaction_features(
        subrecords(mir_records, MIRNA),
        pre_mir_records,
        interaction_list,
        (
            read_2d_structure_loops(MIRNA_LOOPS_FOLDER, precursor_map[mir])
            for mir in interaction_list[MIRNA_FIELD]
        ),
        binding_list,
    )
