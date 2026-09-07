import re


TOTAL_ENERGY_RE = re.compile(r"\((-?\d+(?:\.\d+)?)\)\s*$")
EXTERNAL_LOOP_RE = re.compile(r"^(External loop)\s*:\s*([+-]?\d+)\s*$")
LOOP_RE = re.compile(
    r"^(Interior loop|Hairpin\s+loop|Multi\s+loop)\s*\(\s*(\d+),\s*(\d+)\).*?:\s*([+-]?\d+)\s*$"
)