import os
import re
from collections.abc import Iterable, Iterator
from pathlib import Path

from compgenlnc.consts.fields import PAIR_FIELDS
from compgenlnc.consts.params import LNCRNA, MIRNA
from compgenlnc.consts.paths import BINDING_PREFIX
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import InteractionPair


def read_binding_zone(
    binding_folder: str | os.PathLike, interaction: InteractionPair
) -> BindingZone:
    binding_folder = Path(binding_folder).resolve()
    lnc, mir = interaction[PAIR_FIELDS]
    filename = binding_folder / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
    text = filename.read_text()
    m = re.match(MIRANDA_INFO, text)
    if not m:
        return BindingZone("", "", False)
    lnc_range = range(int(m["lnc_first"]) - 1, int(m["lnc_last"]) - 1)
    mir_range = range(int(m["mir_first"]) - 1, int(m["mir_last"]) - 1)
    return BindingZone(
        lnc, mir, True, float(m["energy"]), lnc_range, mir_range
    )


def iter_binding_zones(
    binding_folder: str | os.PathLike, interactions: Iterable[InteractionPair]
) -> Iterator[BindingZone]:
    return (
        read_binding_zone(binding_folder, interaction)
        for interaction in interactions
    )


def get_subrecord(
    record: SeqSSRecord,
    interaction: InteractionPair,
    binding_zone: BindingZone,
    *,
    rna: str = LNCRNA,
) -> SeqSSRecord:
    if rna == LNCRNA:
        b_range = binding_zone._lnc_range
    elif rna == MIRNA:
        b_range = binding_zone._mir_range
    return SeqSSRecord(
        "_".join(interaction[PAIR_FIELDS]),
        record.seq[b_range.start : b_range.stop],
        record.ss[b_range.start : b_range.stop],
    )


def get_subrecord_list(
    records: Iterable[SeqSSRecord],
    interactions: Iterable[InteractionPair],
    binding_zones: Iterable[BindingZone],
    *,
    rna: str = LNCRNA,
) -> Iterator[SeqSSRecord]:
    return (
        get_subrecord(*args, rna=rna)
        for args in zip(records, interactions, binding_zones)
    )


def get_subrecord_list_from_folder(
    records: Iterable[SeqSSRecord],
    interactions: Iterable[InteractionPair],
    binding_folder: str | os.PathLike,
    *,
    rna: str = LNCRNA,
) -> Iterator[SeqSSRecord]:
    binding_zones = iter_binding_zones(binding_folder, interactions)
    return get_subrecord_list(records, interactions, binding_zones, rna=rna)
