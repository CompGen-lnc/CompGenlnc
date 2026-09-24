import re

EXTERNAL_LOOP_RE = re.compile(
    r"^(?P<kind>External loop)\s*:\s*(?P<energy>[+-]?\d+)\s*$"
)
LOOP_RE = re.compile(
    r"^(?P<kind>Interior loop|Hairpin\s+loop|Multi\s+loop)"
    r"\s*\(\s*(?P<low>\d+),\s*(?P<high>\d+)\).*?:\s*(?P<energy>[+-]?\d+)\s*$"
)
MIRANDA_INFO = re.compile(
    (
        r"Forward:.*?Q:(?P<mir_first>\d+)\s+to\s+(?P<mir_last>\d+).*?"
        r"R:(?P<lnc_first>\d+)\s+to\s+(?P<lnc_last>\d+).*?"
        r"Energy:\s*(?P<energy>-?\d+(?:\.\d+)?)\s*kCal/Mol"
    ),
    re.S,
)
