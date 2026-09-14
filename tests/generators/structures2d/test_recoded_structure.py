import numpy as np
import pytest

from compgenlnc.generators.structures2d import (
    recode_2d_structure,
    recode_2d_structure_list,
)
from compgenlnc.structs import SeqSSRecord, StructureRecord


@pytest.mark.parametrize(
    "ss, expected",
    [
        (
            "((.)..).",
            np.array([6, 6, 5, 7, 5, 5, 7, 5], np.uint8),
        ),
        (
            StructureRecord("..(((.)))"),
            np.array([5, 5, 6, 6, 6, 5, 7, 7, 7], np.uint8),
        ),
        (
            SeqSSRecord("", "", "(((....)))"),
            np.array([6] * 3 + [5] * 4 + [7] * 3, np.uint8),
        ),
    ],
)
def test_recode_2d_structure(ss, expected):
    recoded = recode_2d_structure(ss)
    assert np.array_equal(recoded, expected)


def test_recode_2d_structure_list(ss_fasta_list):
    recoded_list = recode_2d_structure_list(ss_fasta_list)
    for i, (record, expected) in enumerate(zip(ss_fasta_list, recoded_list)):
        assert np.array_equal(recode_2d_structure(record.ss), expected)
