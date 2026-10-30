import compgenlnc as cgl


def test_read_binding_zone(interaction_list, subtests):
    for interaction in interaction_list:
        lnc, mir = interaction[cgl.PAIR_FIELDS]
        folder = cgl.paths.INTERACTIONS_FOLDER / "binding"
        filename = folder / f"{cgl.paths.BINDING_PREFIX}{lnc}_{mir}.dat"
        with subtests.test(i=(lnc, mir)):
            text = filename.read_text()
            binding_zone = cgl.read_binding_zone(folder, interaction)
            assert (text == "No Hits Found above Threshold") == (
                not bool(binding_zone)
            )
