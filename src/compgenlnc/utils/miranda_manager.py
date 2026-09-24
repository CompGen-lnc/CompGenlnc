import os
import re
from pathlib import Path

from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.structs.binding_zone import BindingZone


def get_binding_zone(
    filename: os.PathLike,
) -> BindingZone:
    text = Path(filename).resolve().read_text()
    m = re.match(MIRANDA_INFO, text)
    if not m:
        return BindingZone(False)
    lnc_range = range(int(m["lnc_first"]) - 1, int(m["lnc_last"]) - 1)
    mir_range = range(int(m["mir_first"]) - 1, int(m["mir_last"]) - 1)
    return BindingZone(True, float(m["energy"]), lnc_range, mir_range)