import sys
from pathlib import Path
from typing import Final

# Files
LOG_FILE = sys.stdout

# Seeds
DOC2VEC_MODEL_SEED = 0

# Paths
BASE_PATH: Final[Path] = Path(__file__).resolve().parent.parent
