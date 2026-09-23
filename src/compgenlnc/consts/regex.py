import re

EXTERNAL_LOOP_RE = re.compile(
    r"^(?P<kind>External loop)\s*:\s*(?P<energy>[+-]?\d+)\s*$"
)
LOOP_RE = re.compile(
    r"^(?P<kind>Interior loop|Hairpin\s+loop|Multi\s+loop)"
    r"\s*\(\s*(?P<low>\d+),\s*(?P<high>\d+)\).*?:\s*(?P<energy>[+-]?\d+)\s*$"
)
