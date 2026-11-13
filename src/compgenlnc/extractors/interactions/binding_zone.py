from __future__ import annotations

import multiprocessing as mp
import os
import platform
import re
import subprocess
import tempfile
from collections.abc import Collection, Sequence
from itertools import repeat
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

import compgenlnc.config as config
from compgenlnc.consts.fields import LNCRNA_FIELD, MIRNA_FIELD, PAIR_FIELDS
from compgenlnc.consts.paths import BINDING_PREFIX, MIRANDA_PREFIX
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.dtypes import (
    dtype_interaction_pair,
)
from compgenlnc.fileman.csv_manager import load_interactions
from compgenlnc.fileman.fasta_manager import load_fasta
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import (
    InteractionPair,
)

if TYPE_CHECKING:
    from _typeshed import StrPath


def predict_miranda(
    dir: os.PathLike, lnc_record: SeqSSRecord, mir_record: SeqSSRecord
) -> bool:
    """Predict possible binding zones for an interaction.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory to save the miRanda output file.
    lnc_record : SeqSSRecord
        The record for the lncRNA of the interaction.
    mir_record : SeqSSRecord
        The record for the mature miRNA of the interaction.

    Returns
    -------
    bool
        `True` if the output file was created without interruptions,
        even if there are no binding zones.
    """
    dir = Path(dir).resolve()
    dir.mkdir(parents=True, exist_ok=True)

    filename = dir / f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"

    def to_wsl(path):
        disk, sub_path = str(path).split(":")
        return "/mnt/" + disk.lower() + sub_path.replace("\\", "/")

    def command(mir, lnc, out):
        conda_cmd = ["conda", "run", "-n", config.conda_environment]
        miranda_cmd = ["miranda", mir, lnc, "-out", out]
        return conda_cmd + miranda_cmd

    success = False

    with (
        tempfile.NamedTemporaryFile("w", suffix=".txt") as out_lnc,
        tempfile.NamedTemporaryFile("w", suffix=".txt") as out_mir,
    ):
        out_lnc.write(f">{lnc_record.id}\n")
        out_lnc.write(f"{lnc_record.seq}\n")
        out_lnc.flush()
        out_mir.write(f">{mir_record.id}\n")
        out_mir.write(f"{mir_record.seq}\n")
        out_mir.flush()

        if platform.system() in ["Linux", "Darwin"]:
            args = command(out_mir.name, out_lnc.name, filename)
        elif platform.system() == "Windows":
            cmd = command(
                to_wsl(out_mir.name), to_wsl(out_lnc.name), to_wsl(filename)
            )
            args = ["wsl", "bash", "-i", "-c", " ".join(cmd)]
        else:
            args = []

        try:
            subprocess.run(args, check=True)
            success = True
        except subprocess.CalledProcessError:
            pass

    return success


def _predict_miranda_with_args(args):
    """Wrap `predict_miranda` to use a tuple with arguments."""
    return predict_miranda(*args)


def predict_miranda_from_list(
    dir: os.PathLike,
    lnc_list: Sequence[SeqSSRecord],
    mir_list: Sequence[SeqSSRecord],
    *,
    keep: Collection[InteractionPair] | None = None,
    exclude: Collection[InteractionPair] | None = None,
) -> None:
    """Predict possible binding zones from an interaction set.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory to save the miRanda output file.
    lnc_list : sequence of SeqSSRecord
        Sequence of records for the lncRNAs of the interactions.
    mir_list : sequence of SeqSSRecord
        Sequence of records for the mature miRNAs of the interactions.
    keep : collection of InteractionPair, optional
        Set of interaction pairs to save. Do not apply if `None` or
        empty. Default: `None`.
    exclude : collection of InteractionPair, optional
        Set of interaction pairs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    send_msg("Calculating binding zone for interactions\n")
    lnc_list = list(lnc_list)
    mir_list = list(mir_list)
    args_list = np.fromiter(
        (args for args in zip(repeat(dir), lnc_list, mir_list)),
        dtype=object,
        count=min(len(lnc_list), len(mir_list)),
    )
    mask = np.fromiter(
        (None not in row for row in args_list),
        dtype=np.bool,
        count=len(args_list),
    )
    id_list = np.fromiter(
        ((lnc.id, mir.id) for _, lnc, mir in args_list),
        dtype=dtype_interaction_pair,
        count=len(args_list),
    )
    if keep is not None and len(keep):
        mask &= np.isin(id_list, keep[PAIR_FIELDS])

    total = mask.sum()
    initial = 0

    if exclude is not None and len(exclude):
        mask &= np.isin(id_list, exclude[PAIR_FIELDS], invert=True)
        initial = len(exclude)
    args_list = args_list[mask]

    num_processes = mp.cpu_count() - 2
    with mp.Pool(processes=num_processes) as pool:
        it = pool.imap_unordered(
            _predict_miranda_with_args,
            args_list,
            chunksize=50,
        )

        desc = "  Calculating binding zone"
        it = track_progress(it, desc, total, "pairs", initial)

        for _ in it:
            pass

    send_msg("  Binding zones calculated\n")


def predict_miranda_from_files(
    dir: StrPath,
    interaction_file: StrPath,
    lnc_fasta: StrPath,
    pre_mir_fasta: StrPath,
    mat_mir_fasta: StrPath,
    *,
    keep: Collection[InteractionPair] | None = None,
    exclude: Collection[InteractionPair] | None = None,
) -> None:
    """Predict possible binding zones from the source files.

    Parameters
    ----------
    dir : str or os.PathLike
        Path of the directory to save the miRanda output file.
    interaction_file : str or os.PathLike
        Path of the CSV file with the interaction pairs.
    lnc_fasta : str or os.PathLike
        Path of the FASTA file with the lncRNA sequences.
    pre_mir_fasta : str or os.PathLike
        Path of the FASTA file with the precursor miRNA sequences.
    mat_mir_fasta : str or os.PathLike
        Path of the FASTA file with the mature miRNA sequences.
    keep : collection of InteractionPair, optional
        Set of interaction pairs to save. Do not apply if `None` or
        empty. Default: `None`.
    exclude : collection of InteractionPair, optional
        Set of interaction pairs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    interactions = load_interactions(interaction_file)

    lnc_records = load_fasta(lnc_fasta)
    mir_records = load_fasta(pre_mir_fasta) + load_fasta(mat_mir_fasta)
    lnc_map = {rec.id: rec for rec in lnc_records}
    mir_map = {rec.id: rec for rec in mir_records}
    lnc_list = np.select(
        [interactions[LNCRNA_FIELD] == key for key in lnc_map.keys()],
        lnc_map.values(),
        None,
    )
    mir_list = np.select(
        [interactions[MIRNA_FIELD] == key for key in mir_map.keys()],
        mir_map.values(),
        None,
    )
    predict_miranda_from_list(
        dir, lnc_list, mir_list, keep=keep, exclude=exclude
    )


def extract_binding_zone(
    miranda_dir: StrPath,
    lnc: str,
    mir: str,
    binding_dir: StrPath,
) -> None:
    """Extract the lowest energy binding zone for an interaction.

    Parameters
    ----------
    miranda_dir : str or os.PathLike
        Path of the directory containing the miRanda output file.
    lnc : str
        ID of the lncRNA.
    mir : str
        ID of the mature miRNA.
    binding_dir : str or os.PathLike
        Path of the directory where the binding zone will be saved.
    """
    miranda_dir = Path(miranda_dir).resolve()
    binding_dir = Path(binding_dir).resolve()

    binding_dir.mkdir(parents=True, exist_ok=True)

    miranda_file = miranda_dir / f"{MIRANDA_PREFIX}{lnc}_{mir}.dat"
    text = miranda_file.read_text()
    hits = list(re.finditer(MIRANDA_INFO, text))
    if hits:
        min_energy_hit = min(hits, key=lambda hit: hit["energy"])[0]
    else:
        min_energy_hit = "No Hits Found above Threshold"
    filename = binding_dir / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
    with open(filename, "w") as out_file:
        out_file.write(min_energy_hit)


def extract_binding_zone_from_dir(
    miranda_dir: StrPath,
    binding_dir: StrPath,
    *,
    lnc_set: Collection[str] | None = None,
    mir_set: Collection[str] | None = None,
    keep: Collection[InteractionPair] | None = [],
    exclude: Collection[InteractionPair] | None = [],
) -> None:
    """Extract the lowest energy binding zone from miRanda output files.

    Parameters
    ----------
    miranda_dir : str or os.PathLike
        Path of the directory containing the miRanda output files.
    binding_dir : str or os.PathLike
        Path of the directory where the binding zones will be saved.
    lnc_set : collection of str, optional
        Set of IDs of the lncRNA in the interactions.
    mir_set : collection of str, optional
        Set of IDs of the mature miRNA in the interactions.
    keep : collection of InteractionPair, optional
        Set of interaction pairs to save. Do not apply if `None` or
        empty. Default: `None`.
    exclude : collection of InteractionPair, optional
        Set of interaction pairs to ignore. Do not apply if `None` or
        empty. Default: `None`.
    """
    miranda_files = os.listdir(miranda_dir)
    keep = ["_".join(pair[PAIR_FIELDS]) for pair in keep]
    exclude = ["_".join(pair[PAIR_FIELDS]) for pair in exclude]
    for file in miranda_files:
        prefix = file.find(MIRANDA_PREFIX)
        ext = file.find(".dat")
        if prefix < 0 or ext < 0:
            continue

        interaction = file[len(MIRANDA_PREFIX) : ext]
        if (
            keep
            and interaction not in keep
            or exclude
            and interaction in exclude
        ):
            continue
        [lnc, mir] = interaction.split("_")
        if lnc_set and lnc not in lnc_set or mir_set and mir not in mir_set:
            continue
        extract_binding_zone(miranda_dir, lnc, mir, binding_dir)
