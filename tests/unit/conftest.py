from pathlib import Path

import compgenlnc as cgl

cgl.set_base_path(Path(__file__).resolve().parent)

pytest_plugins = [
    "fixtures.deterministic_models",
    "fixtures.id_set",
    "fixtures.interaction_set",
    "fixtures.seq_ss_set",
    "fixtures.sequence_set",
    "fixtures.structure_set",
    "fixtures.index_set",
]
