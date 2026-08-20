import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import *
from compgenlnc.extractors.mirna.fasta_constructor import iter_mirna


def join_seq_ss(seq: str, ss: str, mir: str) -> None:
    path = MIRNA_SEQ_SS_FOLDER / f'{SEQ_SS_PREFIX}{mir}.txt'
    with open(path, 'w') as out_file:
        out_file.write(seq + '\n')
        out_file.write(ss)


def fasta_to_seq_ss() -> int:
    mir_dict_seq = {
        record.id: record.seq.replace('U', 'T')
        for record in iter_mirna(PRE_MIRNA_FASTA)
    }
    mir_list = list(mir_dict_seq.keys())
    MIRNA_SEQ_SS_FOLDER.mkdir(parents=True, exist_ok=True)
    total = 0

    for mir in mir_list:
        try:
            seq = mir_dict_seq[mir]
            join_seq_ss(seq, '', mir)
            total += 1
        except:
            print(mir)

    return total
