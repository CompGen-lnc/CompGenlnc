import compgenlnc as cgl


def test_read_2d_structure_loops(seq_ss_list):
    for record in seq_ss_list:
        loop_counter = cgl.read_2d_structure_loops(
            cgl.paths.MIRNA_LOOPS_FOLDER, record.id
        )
        total_energy = loop_counter.loops[cgl.ENERGY_FIELD].sum().item() / 100
        assert loop_counter.energy == total_energy


def test_read_2d_structure_loops_from_folder(example_id_filter):
    loops_iter = cgl.read_2d_structure_loops_from_folder(
        cgl.paths.MIRNA_LOOPS_FOLDER, keep=example_id_filter
    )
    for id_, loop_counter in zip(example_id_filter, loops_iter):
        expected = cgl.read_2d_structure_loops(
            cgl.paths.MIRNA_LOOPS_FOLDER, id_
        )
        assert loop_counter == expected
