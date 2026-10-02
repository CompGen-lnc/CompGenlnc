from compgenlnc.consts import PAIR_FIELDS
from compgenlnc.config.paths import BINDING_PREFIX
from compgenlnc.fileman import get_binding_zone

from constants import INTERACTIONS_FOLDER


def test_get_binding_zone(interaction_list, subtests):
    for interaction in interaction_list:
        lnc, mir = interaction[PAIR_FIELDS]
        folder = INTERACTIONS_FOLDER / "binding"
        filename = folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        with subtests.test(i=(lnc, mir)):
            text = filename.read_text()
            binding_zone = get_binding_zone(folder, interaction)
            assert (text == "No Hits Found above Threshold") == (
                not bool(binding_zone)
            )
