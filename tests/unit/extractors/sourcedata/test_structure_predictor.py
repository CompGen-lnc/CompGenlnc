import pytest

import compgenlnc as cgl


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            (
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA"
            ),
            cgl.structs.StructureRecord(
                "(((((.(((((((((((((((((((((.....(((...(("
                "((....)))).)))))))))))))))))))))))))))))"
            ),
        ),
        (
            cgl.structs.SeqRecord(
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU"
            ),
            cgl.structs.StructureRecord(
                "(((..(((.(((.(((((((((((((.........("
                "((......)))))))))))))))).))).))).)))"
            ),
        ),
        (
            (
                "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCC"
                "CUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU"
            ),
            cgl.structs.StructureRecord(
                "(((.(((((((((((((((((((((((((((...)))"
                "))).........))))))))))))))))))))).)))"
            ),
        ),
    ],
)
def test_extract_2d_structure(seq, expected):
    ss = cgl.extract_2d_structure(seq)
    assert len(seq) == len(ss)
    assert ss == expected


@pytest.mark.parametrize(
    "seq, expected",
    [
        (
            cgl.structs.SeqSSRecord(
                "hsa-let-7a-1",
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
                "",
            ),
            cgl.structs.SeqSSRecord(
                "hsa-let-7a-1",
                "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
                "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
                "(((((.(((((((((((((((((((((.....(((...(("
                "((....)))).)))))))))))))))))))))))))))))",
            ),
        ),
        (
            cgl.structs.SeqSSRecord(
                "hsa-let-7a-2",
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU",
                "",
            ),
            cgl.structs.SeqSSRecord(
                "hsa-let-7a-2",
                "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACA"
                "UCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU",
                "(((..(((.(((.(((((((((((((.........("
                "((......)))))))))))))))).))).))).)))",
            ),
        ),
        (
            cgl.structs.SeqSSRecord(
                "hsa-let-7a-3",
                "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCC"
                "CUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU",
                "",
            ),
            cgl.structs.SeqSSRecord(
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
    ss = cgl.extract_2d_structure_identified(seq)
    assert len(seq.seq) == len(ss.ss)
    assert ss == expected


def test_extract_2d_structure_list(seq_ss_list):
    ss_list = cgl.extract_2d_structure_list(seq_ss_list)
    for i, expected in enumerate(ss_list):
        ss = cgl.extract_2d_structure(seq_ss_list[i])
        assert ss == expected
