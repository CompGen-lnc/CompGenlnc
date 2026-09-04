import os
from pathlib import Path
from typing import Iterable, Iterator

from compgenlnc.config.paths import (
    SEQ_SS_PREFIX,
)
from compgenlnc.structs import SeqSSRecord


def load_fasta(
        filename: str | os.PathLike,
        /, *,
        mode: str = 'seq',
        keep: Iterable[str] | None = None,
) -> Iterator[SeqSSRecord]:
    """Load a FASTA file and iterate the RNA molecules into them.

    Args:
        filename: Path of the FASTA file.
        mode: Type of the data saved in the file. If 'seq', it will treat the
            data as sequences. If 'ss', it will treat the data as 2D
            structures.

    Yields:
        SeqSSRecord: A record with the sequence or structure and its name.

    Raises:
        ValueError: If mode is any other than 'seq' or 'ss'.
    """
    if mode not in ['seq', 'ss']:
        raise ValueError(f"Value {mode} for mode is not allowed")
    
    filename = Path(filename).resolve()
    id_ = seq = ss = value = ''
    with open(filename) as file:
        # Add the '>' at the end for capturing the last molecule
        for raw_line in file.readlines() + ['>']:
            if raw_line[0] != '>':
                if id_ and raw_line.strip():
                    value += raw_line.split()[0]
                continue

            if mode == 'seq':
                seq = value
            elif mode == 'ss':
                ss = value

            if id_:    
                yield SeqSSRecord(id_, seq, ss)
                id_ = ''
                value = ''

            if not raw_line[1:].strip():
                continue

            line = raw_line[1:].split()[0]
            if keep is None or line in keep:
                id_ = line



def save_fasta(
        filename: str | os.PathLike,
        molecules: Iterable[SeqSSRecord],
        /, *,
        mode: str,
) -> None:
    """Save sequences or 2D structures from a group of RNA molecules into a
    FATSA file.

    Args:
        filename: Path of the FASTA file.
        molecules: Group of RNA molecules to 
        mode: Type of the data saved in the file. If 'seq', it will save the
            sequences. If 'ss', it will save the 2Dstructures.

    Raises:
        ValueError: If mode is any other than 'seq' or 'ss'.
    """
    if mode not in ['seq', 'ss']:
        raise ValueError(f"Value {mode} for mode is not allowed")
    
    filename = Path(filename).resolve()
    filename.parent.mkdir(parents=True, exist_ok=True)
    with open(filename, 'w') as out_file:
        for record in molecules:
            value = record.seq if mode == 'seq' else record.ss
            out_file.write(f'>{record.id}\n')
            if not value:
                continue
            out_file.write(f'{'\n'.join([
                value[i : i + 60] for i in range(0, len(value), 60)
            ])}\n')

    
def iter_seq_fasta(
        filename: str | os.PathLike,
        filter: Iterable[str] | None = None
) -> Iterator[SeqSSRecord]:
    """Iterate all the RNA sequences from a FASTA file that fit the 
    filter.

    Args:
        filename: Path of the FASTA file.
        filter: A list with all substrings contained in the name of the 
            desired sequences. If None, it is the same as using `load_fasta`
            with mode as 'seq'.

    Yields:
        SeqSSRecord: A record with the sequence and its name.
    """
    filename = Path(filename).resolve()

    for record in load_fasta(filename, mode='seq'):
        if (filter and
            not any(prefix in record.id for prefix in filter)):
            continue
        yield record


def read_sequence(
        filename: str | os.PathLike,
        mature_only: bool = False
) -> str:
    """Read the first line of a seq+ss file and normalize to ACGT
    alphabet.

    Args:
        filename: Path of the seq+ss file.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    filename = Path(filename).resolve()

    with open(filename) as f:
        seq = f.readlines()[0].strip()
    if mature_only:
        seq = ''.join(c for c in seq if c.isupper())
    return seq.upper().replace('U', 'T')


def read_structure(
        filename: str | os.PathLike,
        mature_only: bool = False
) -> str:
    """Read the first line of a seq+ss file and normalize to ACGT
    alphabet.

    Args:
        filename: Path of the seq+ss file.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    filename = Path(filename).resolve()

    with open(filename) as f:
        ss = f.readlines()[1].strip()
    if mature_only:
        ss = ''.join(c for c in ss if c.isupper())
    return ss.upper().replace('U', 'T')


def iter_seq_ss(
        folder: str | os.PathLike,
        keep: Iterable[str] | None = None,
        mature_only: bool = False
) -> Iterator[SeqSSRecord]:
    """Iterate all the seq+ss files from the same folder that fit the
    filter.

    Args:
        folder: Path of the folder containing seq+ss files.
        keep: List of the RNA names to iterate. Iterate all if None.
        mature_only: Only saves the mature section of the sequence if
            True.
    """
    folder = Path(folder).resolve()

    for filename in sorted(os.listdir(folder)):
        if not filename.startswith(SEQ_SS_PREFIX):
            continue
        id_ = os.path.splitext(filename[len(SEQ_SS_PREFIX):])[0]
        if keep and id_ not in keep:
            continue
        seq = read_sequence(folder / filename, mature_only)
        ss = read_structure(folder / filename, mature_only)
        yield SeqSSRecord(id_, seq, ss)