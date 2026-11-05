import numpy as np

import compgenlnc as cgl


def test_get_substructure_loops(
    interaction_list, pre_mirna_id_list, precursors_consverter, subtests
):
    interactions = interaction_list[
        ~np.isin(interaction_list[cgl.MIRNA_FIELD], pre_mirna_id_list)
    ]
    for interaction in interactions[cgl.PAIR_FIELDS]:
        lnc, mir = interaction
        precursor = precursors_consverter[mir]
        binding_folder = cgl.paths.INTERACTIONS_FOLDER / "binding"
        with subtests.test(i=(lnc, mir)):
            mir_loops = cgl.read_2d_structure_loops(
                cgl.paths.MIRNA_LOOPS_FOLDER, precursor
            )
            binding_zone = cgl.read_binding_zone(binding_folder, interaction)
            sub_loops = cgl.get_substructure_loops(
                mir_loops, binding_zone, molecule=cgl.MIRNA
            )
            start, end = binding_zone.mir_start, binding_zone.mir_end
            mask = (mir_loops.loops[cgl.LOW_POSITION] >= start) & (
                mir_loops.loops[cgl.HIGH_POSITION] <= end
            )
            expected = mir_loops.loops[mask]
            assert np.array_equal(sub_loops.loops, expected)


def test_get_substructure_loops_from_list(
    interaction_list, pre_mirna_id_list, precursors_consverter
):
    loops_list = []
    binding_list = []
    interactions = interaction_list[
        ~np.isin(interaction_list[cgl.MIRNA_FIELD], pre_mirna_id_list)
    ]
    for interaction in interactions[cgl.PAIR_FIELDS]:
        _, mir = interaction
        precursor = precursors_consverter[mir]
        binding_folder = cgl.paths.BINDING_FOLDER
        mir_loops = cgl.read_2d_structure_loops(
            cgl.paths.MIRNA_LOOPS_FOLDER, precursor
        )
        binding_zone = cgl.read_binding_zone(binding_folder, interaction)
        loops_list.append(mir_loops)
        binding_list.append(binding_zone)

    new_loops_list = cgl.get_substructure_loops_from_list(
        loops_list, binding_list, molecule=cgl.MIRNA
    )
    for sub_loops, mir_loops, binding_zone in zip(
        new_loops_list, loops_list, binding_list
    ):
        expected = cgl.get_substructure_loops(
            mir_loops, binding_zone, molecule=cgl.MIRNA
        )
        assert sub_loops == expected


def test_gen_substructure_loops_dict(
    interaction_list, pre_mirna_id_list, precursors_consverter, tmp_path
):
    dict_path = tmp_path / "sub_loops.dict"
    loops_list = []
    binding_list = []
    interactions = interaction_list[
        ~np.isin(interaction_list[cgl.MIRNA_FIELD], pre_mirna_id_list)
    ]
    for interaction in interactions[cgl.PAIR_FIELDS]:
        _, mir = interaction
        precursor = precursors_consverter[mir]
        binding_folder = cgl.paths.BINDING_FOLDER
        mir_loops = cgl.read_2d_structure_loops(
            cgl.paths.MIRNA_LOOPS_FOLDER, precursor
        )
        binding_zone = cgl.read_binding_zone(binding_folder, interaction)
        loops_list.append(mir_loops)
        binding_list.append(binding_zone)

    cgl.gen_substructure_loops_dict(
        dict_path, interaction_list, loops_list, binding_list, molecule=cgl.MIRNA
    )
    loops_dict = cgl.load_dict(dict_path)
    for interaction, mir_loops, binding_zone in zip(
        interactions, loops_list, binding_list
    ):
        if interaction[cgl.IS_POSITIVE_FIELD]:
            expected = cgl.get_substructure_loops(
                mir_loops, binding_zone, molecule=cgl.MIRNA
            ).loops_vectorized
        else:
            expected = np.array([-1] * 3 + [0] * 4)
        vector = loops_dict["_".join(interaction[cgl.PAIR_FIELDS])]
        assert np.array_equal(vector, expected)
