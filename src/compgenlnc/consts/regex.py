import re


TOTAL_ENERGY_RE = re.compile(r"\((-?\d+(?:\.\d+)?)\)\s*$")
LOOP_RE = re.compile(
    r"^(External loop|Interior loop|Hairpin\s+loop|Multi\s+loop)\s*\(\s*(\d+),\s*(\d+)\).*?:\s*([+-]?\d+)\s*$"
)