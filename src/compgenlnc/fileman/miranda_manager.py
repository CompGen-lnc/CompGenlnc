import os
import re
from pathlib import Path

from compgenlnc.consts.fields import PAIR_FIELDS
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.config.paths import BINDING_PREFIX
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.typing.numpy_dtypes import InteractionPair


def get_binding_zone(
    folder: str | os.PathLike, interaction: InteractionPair
) -> BindingZone:
    folder = Path(folder).resolve()
    lnc, mir = interaction[PAIR_FIELDS]
    filename = folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
    text = filename.read_text()
    m = re.match(MIRANDA_INFO, text)
    if not m:
        return BindingZone(False)
    lnc_range = range(int(m["lnc_first"]) - 1, int(m["lnc_last"]) - 1)
    mir_range = range(int(m["mir_first"]) - 1, int(m["mir_last"]) - 1)
    return BindingZone(True, float(m["energy"]), lnc_range, mir_range)
