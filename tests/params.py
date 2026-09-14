from compgenlnc.structs import SeqSSRecord

id_list = [
    "hsa-let-7a-1",
    "hsa-let-7a-2",
    "hsa-let-7a-3",
    "hsa-let-7b",
    "hsa-let-7c",
    "hsa-let-7d",
    "hsa-let-7e",
    "hsa-let-7f-1",
    "hsa-let-7f-2",
    "hsa-let-7g",
]

seq_list = [
    SeqSSRecord(
        "hsa-let-7a-1",
        (
            "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACC"
            "CACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA"
        ),
        "",
    ),
    SeqSSRecord(
        "hsa-let-7a-2",
        (
            "AGGUUGAGGUAGUAGGUUGUAUAGUUUAGAAUUACAUCAAGGGAGAUAACUGUACAGCCUCCUAGCUUUCCU"
        ),
        "",
    ),
    SeqSSRecord(
        "hsa-let-7a-3",
        (
            "GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCCCUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU"
        ),
        "",
    ),
    SeqSSRecord(
        "hsa-let-7b",
        (
            "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUG"
            "UUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"
        ),
        "",
    ),
    SeqSSRecord(
        "hsa-let-7c",
        (
            "GCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACA"
            "CCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC"
        ),
        "",
    ),
]

T_converted_seq_list = [
    SeqSSRecord(seq.id, seq.seq.convert_to_T(), seq.ss) for seq in seq_list
]
