import pytest

from compgenlnc.extractors.structures2d import (
    extract_2d_structure,
    extract_2d_structure_identified,
    extract_2d_structure_list,
)
from compgenlnc.structs import SeqRecord, SeqSSRecord, StructureRecord


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            (
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA"
            ),
            StructureRecord(
                "(((((.(((((((((((((((((((((.....(((...(("
                "((....)))).)))))))))))))))))))))))))))))"
            ),
        ),
        (
            SeqRecord(
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU"
            ),
            StructureRecord(
                "(((..(((.(((.(((((((((((((.........("
                "((......)))))))))))))))).))).))).)))"
            ),
        ),
        (
            (
                "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCC"
                "CUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU"
            ),
            StructureRecord(
                "(((.(((((((((((((((((((((((((((...)))"
                "))).........))))))))))))))))))))).)))"
            ),
        ),
    ],
)
def test_extract_2d_structure(seq, expected):
    ss = extract_2d_structure(seq)
    assert len(seq) == len(ss)
    assert ss == expected


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            SeqSSRecord(
                "hsa-let-7a-1",
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
                "",
            ),
            SeqSSRecord(
                "hsa-let-7a-1",
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
                "(((((.(((((((((((((((((((((.....(((...(("
                "((....)))).)))))))))))))))))))))))))))))",
            ),
        ),
        (
            SeqSSRecord(
                "hsa-let-7a-2",
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU",
                "",
            ),
            SeqSSRecord(
                "hsa-let-7a-2",
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU",
                "(((..(((.(((.(((((((((((((.........("
                "((......)))))))))))))))).))).))).)))",
            ),
        ),
        (
            SeqSSRecord(
                "hsa-let-7a-3",
                "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCC"
                "CUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU",
                "",
            ),
            SeqSSRecord(
                "hsa-let-7a-3",
                "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCC"
                "CUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU",
                "(((.(((((((((((((((((((((((((((...)))"
                "))).........))))))))))))))))))))).)))",
            ),
        ),
    ],
)
def test_extract_2d_structure_identified(seq, expected):
    ss = extract_2d_structure_identified(seq)
    assert len(seq.seq) == len(ss.ss)
    assert ss == expected


def test_extract_2d_structure_list(seq_ss_filtered_list):
    ss_list = extract_2d_structure_list(seq_ss_filtered_list)
    for i, expected in enumerate(ss_list):
        ss = extract_2d_structure(seq_ss_filtered_list[i])
        assert ss == expected
