import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord

from compgenlnc.config.paths import (
    LNCRNA_FASTA,
    LNCRNA_SEQ_SS_FOLDER,
    PRE_MIRNA_FASTA,
    MIRNA_SEQ_SS_FOLDER,
    SEQ_SS_PREFIX,
)
from compgenlnc.utils.fasta_manager import iter_seq_fasta


def join_seq_ss(
        seq: str,
        ss: str,
        id_: str,
        folder: str | os.PathLike
) -> None:
    folder = Path(folder).resolve()

    path = folder / f'{SEQ_SS_PREFIX}{id_}.dat'
    with open(path, 'w') as out_file:
        out_file.write(seq + '\n')
        out_file.write(ss)


def fasta_to_seq_ss(
        seq_file: str | os.PathLike,
        ss_file: str | os.PathLike,
        folder: str | os.PathLike
) -> int:
    seq_file = Path(seq_file).resolve()
    folder = Path(folder).resolve()

    seq_dict = {
        record.id: str(record.seq).replace('U', 'T')
        for record in iter_seq_fasta(seq_file)
    }
    id_list = list(seq_dict.keys())
    id_list.sort()
    folder.mkdir(parents=True, exist_ok=True)
    total = 0

    for id_ in id_list:
        try:
            seq = seq_dict[id_]
            join_seq_ss(seq, '', id_, folder)
            total += 1
        except:
            print(id_)

    return total


def mirna_to_seq_ss() -> int:
    return fasta_to_seq_ss(PRE_MIRNA_FASTA, MIRNA_SEQ_SS_FOLDER)


def lncrna_to_seq_ss() -> int:
    return fasta_to_seq_ss(LNCRNA_FASTA, LNCRNA_SEQ_SS_FOLDER)
