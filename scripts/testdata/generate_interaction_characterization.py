import os
import sys
from pathlib import Path

from compgenlnc.generators.interactions import (
    predict_miranda_from_files,
    extract_binding_zone_from_folder,
)

raiz_proyecto = str(Path(__file__).resolve().parents[2])
if raiz_proyecto not in sys.path:
    sys.path.append(raiz_proyecto)

from tests.constants import (
    INTERACTIONS_CSV,
    INTERACTIONS_FOLDER,
    LNCRNA_FASTA,
    MAT_MIRNA_FASTA,
    PRE_MIRNA_FASTA,
)


def generate_binding_files():
    miranda_folder = INTERACTIONS_FOLDER / "miranda"
    binding_folder = INTERACTIONS_FOLDER / "binding"
    if not miranda_folder.exists():
        predict_miranda_from_files(
            INTERACTIONS_CSV,
            LNCRNA_FASTA,
            PRE_MIRNA_FASTA,
            MAT_MIRNA_FASTA,
            miranda_folder,
        )
    if not binding_folder.exists():
        extract_binding_zone_from_folder(miranda_folder, binding_folder)


if __name__ == "__main__":
    generate_binding_files()
