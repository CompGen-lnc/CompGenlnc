import numpy as np
import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "ss, expected",
    [
        (
            "((.)..).",
            np.array([6, 6, 5, 7, 5, 5, 7, 5], np.uint8),
        ),
        (
            cgl.structs.StructureRecord("..(((.)))"),
            np.array([5, 5, 6, 6, 6, 5, 7, 7, 7], np.uint8),
        ),
        (
            cgl.structs.SeqSSRecord("", "", "(((....)))"),
            np.array([6] * 3 + [5] * 4 + [7] * 3, np.uint8),
        ),
    ],
)
def test_recode_2d_structure(ss, expected):
    recoded = cgl.recode_2d_structure(ss)
    assert np.array_equal(recoded, expected)


def test_recode_2d_structure_list(ss_fasta_list):
    recoded_list = cgl.recode_2d_structure_list(ss_fasta_list)
    for i, (record, expected) in enumerate(zip(ss_fasta_list, recoded_list)):
        assert np.array_equal(cgl.recode_2d_structure(record.ss), expected)
