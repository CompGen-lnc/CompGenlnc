from pathlib import Path

import compgenlnc as cgl

if __name__ == "__main__":
    unit_tests_path = "tests/unit"
    cgl.set_base_path(unit_tests_path)
    basefiles = "basefiles/test"
    cgl.load_source_data(basefiles)
    cgl.gen_source_data()
    cgl.gen_database()
