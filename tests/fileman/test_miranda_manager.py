from compgenlnc.consts import PAIR_FIELDS
from compgenlnc.config.paths import BINDING_PREFIX
from compgenlnc.fileman import get_binding_zone

from constants import INTERACTIONS_FOLDER


def test_get_binding_zone(interaction_list, subtests):
    for lnc, mir in interaction_list[PAIR_FIELDS]:
        filename = (
            INTERACTIONS_FOLDER
            / "binding"
            / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
        )
        with subtests.test(i=(lnc, mir)):
            text = filename.read_text()
            binding_zone = get_binding_zone(filename)
            assert (text == "No Hits Found above Threshold") == (
                not bool(binding_zone)
            )
