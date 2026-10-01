import numpy as np

from compgenlnc.config.paths import BINDING_PREFIX, LOOPS_PREFIX
from compgenlnc.consts.params import MIRNA_TYPE
from compgenlnc.extractors.interactions import (
    gen_substructure_loops_dict,
    get_substructure_loops,
    get_substructure_loops_from_list,
)
from compgenlnc.fileman import (
    get_2d_struture_loops,
    get_binding_zone,
    load_dict,
)

from constants import INTERACTIONS_FOLDER, SS_LOOPS_FOLDER


def test_get_substructure_loops(
    interaction_list, pre_mirna_id_list, precursors_consverter, subtests
):
    interactions = interaction_list[
        ~np.isin(interaction_list["miRNA"], pre_mirna_id_list)
    ]
    for lnc, mir in interactions[["lncRNA", "miRNA"]]:
        precursor = precursors_consverter[mir]
        mir_loops_path = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{precursor}.dat"
        binding_zone_path = (
            INTERACTIONS_FOLDER
            / "binding"
            / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        )
        with subtests.test(i=(lnc, mir)):
            mir_loops = get_2d_struture_loops(mir_loops_path)
            binding_zone = get_binding_zone(binding_zone_path)
            sub_loops = get_substructure_loops(
                mir_loops, binding_zone, MIRNA_TYPE
            )
            start, end = binding_zone.mir_start, binding_zone.mir_end
            mask = (mir_loops.loops["low"] >= start) & (
                mir_loops.loops["high"] <= end
            )
            expected = mir_loops.loops[mask]
            assert np.array_equal(sub_loops.loops, expected)


def test_get_substructure_loops_from_list(
    interaction_list, pre_mirna_id_list, precursors_consverter
):
    loops_list = []
    binding_list = []
    interactions = interaction_list[
        ~np.isin(interaction_list["miRNA"], pre_mirna_id_list)
    ]
    for lnc, mir in interactions[["lncRNA", "miRNA"]]:
        precursor = precursors_consverter[mir]
        mir_loops_path = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{precursor}.dat"
        binding_zone_path = (
            INTERACTIONS_FOLDER
            / "binding"
            / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        )
        mir_loops = get_2d_struture_loops(mir_loops_path)
        binding_zone = get_binding_zone(binding_zone_path)
        loops_list.append(mir_loops)
        binding_list.append(binding_zone)

    new_loops_list = get_substructure_loops_from_list(
        loops_list, binding_list, MIRNA_TYPE
    )
    for sub_loops, mir_loops, binding_zone in zip(
        new_loops_list, loops_list, binding_list
    ):
        expected = get_substructure_loops(mir_loops, binding_zone, MIRNA_TYPE)
        assert sub_loops == expected


def test_gen_substructure_loops_dict(
    interaction_list, pre_mirna_id_list, precursors_consverter, tmp_path
):
    dict_path = tmp_path / "sub_loops.dict"
    loops_list = []
    binding_list = []
    interactions = interaction_list[
        ~np.isin(interaction_list["miRNA"], pre_mirna_id_list)
    ]
    for lnc, mir in interactions[["lncRNA", "miRNA"]]:
        precursor = precursors_consverter[mir]
        mir_loops_path = SS_LOOPS_FOLDER / f"{LOOPS_PREFIX}{precursor}.dat"
        binding_zone_path = (
            INTERACTIONS_FOLDER
            / "binding"
            / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        )
        mir_loops = get_2d_struture_loops(mir_loops_path)
        binding_zone = get_binding_zone(binding_zone_path)
        loops_list.append(mir_loops)
        binding_list.append(binding_zone)

    gen_substructure_loops_dict(
        dict_path, interaction_list, loops_list, binding_list, MIRNA_TYPE
    )
    loops_dict = load_dict(dict_path)
    for vector, mir_loops, binding_zone in zip(
        loops_dict.values(), loops_list, binding_list
    ):
        expected = get_substructure_loops(mir_loops, binding_zone, MIRNA_TYPE)
        assert np.array_equal(vector, expected.loops_vectorized)
