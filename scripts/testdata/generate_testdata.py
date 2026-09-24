from generate_interaction_characterization import generate_binding_files
from generate_sequence_characterization import (
    generate_lnc_seq_characterization,
    generate_mat_mir_seq_characterization,
    generate_pre_mir_seq_characterization,
)


if __name__ == "__main__":
    generate_lnc_seq_characterization()
    generate_pre_mir_seq_characterization()
    generate_mat_mir_seq_characterization()

    generate_binding_files()
