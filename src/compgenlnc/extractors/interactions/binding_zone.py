import multiprocessing as mp
import os
import platform
import re
import subprocess
import tempfile
from collections.abc import Iterable
from itertools import repeat
from pathlib import Path

import numpy as np

from compgenlnc.consts.fields import LNCRNA_FIELD, MIRNA_FIELD, PAIR_FIELDS
from compgenlnc.consts.paths import BINDING_PREFIX, MIRANDA_PREFIX
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.fileman.csv_manager import load_interactions
from compgenlnc.fileman.fasta_manager import load_fasta
from compgenlnc.fileman.log_manager import send_msg, track_progress
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import (
    dtype_interaction_pair,
    InteractionPair,
)


def predict_miranda(
    lnc_record: SeqSSRecord, mir_record: SeqSSRecord, folder: os.PathLike
) -> bool:
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)

    filename = folder / f"{MIRANDA_PREFIX}{lnc_record.id}_{mir_record.id}.dat"

    def to_wsl(path):
        disk, sub_path = str(path).split(":")
        return "/mnt/" + disk.lower() + sub_path.replace("\\", "/")

    def command(mir, lnc, out):
        return ["conda", "run", "-n", "base", "miranda", mir, lnc, "-out", out]

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
        except subprocess.CalledProcessError as e:
            pass

    return success


def _predict_miranda_with_args(args):
    return predict_miranda(*args)


def predict_miranda_from_list(
    lnc_list: Iterable[SeqSSRecord],
    mir_list: Iterable[SeqSSRecord],
    folder: os.PathLike,
    *,
    keep: Iterable[InteractionPair] | None = None,
    exclude: Iterable[InteractionPair] | None = None,
) -> None:
    send_msg("Calculating binding zone for interactions\n")
    lnc_list = list(lnc_list)
    mir_list = list(mir_list)
    args_list = np.fromiter(
        (args for args in zip(lnc_list, mir_list, repeat(folder))),
        dtype=object,
        count=min(len(lnc_list), len(mir_list)),
    )
    mask = np.fromiter(
        (None not in row for row in args_list),
        dtype=np.bool,
        count=len(args_list),
    )
    id_list = np.fromiter(
        ((lnc.id, mir.id) for lnc, mir, _ in args_list),
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
        it = track_progress(it, desc=desc, total=total, initial=initial)

        for _ in it:
            pass

    send_msg("  Binding zones calculated\n")


def predict_miranda_from_files(
    interaction_file: str | os.PathLike,
    lnc_fasta: str | os.PathLike,
    pre_mir_fasta: str | os.PathLike,
    mat_mir_fasta: str | os.PathLike,
    folder: str | os.PathLike,
    *,
    keep: Iterable[InteractionPair] | None = None,
    exclude: Iterable[InteractionPair] | None = None,
) -> None:
    interactions = load_interactions(interaction_file)

    lnc_records = list(load_fasta(lnc_fasta))
    mir_records = list(load_fasta(pre_mir_fasta)) + list(
        load_fasta(mat_mir_fasta)
    )
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
        lnc_list, mir_list, folder, keep=keep, exclude=exclude
    )


def extract_binding_zone(
    lnc: str,
    mir: str,
    miranda_folder: os.PathLike,
    binding_folder: os.PathLike,
) -> None:
    miranda_folder = Path(miranda_folder).resolve()
    binding_folder = Path(binding_folder).resolve()

    binding_folder.mkdir(parents=True, exist_ok=True)

    miranda_file = miranda_folder / f"{MIRANDA_PREFIX}{lnc}_{mir}.dat"
    text = miranda_file.read_text()
    hits = list(re.finditer(MIRANDA_INFO, text))
    if hits:
        min_energy_hit = min(hits, key=lambda hit: hit["energy"])[0]
    else:
        min_energy_hit = "No Hits Found above Threshold"
    filename = binding_folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
    with open(filename, "w") as out_file:
        out_file.write(min_energy_hit)


def extract_binding_zone_from_folder(
    miranda_folder: os.PathLike,
    binding_folder: os.PathLike,
    /,
    lnc_list: Iterable[str] | None = None,
    mir_list: Iterable[str] | None = None,
    *,
    keep: Iterable[InteractionPair] | None = [],
    exclude: Iterable[InteractionPair] | None = [],
) -> None:
    miranda_files = os.listdir(miranda_folder)
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
        if (
            lnc_list
            and lnc not in lnc_list
            or mir_list
            and mir not in mir_list
        ):
            continue
        extract_binding_zone(lnc, mir, miranda_folder, binding_folder)
