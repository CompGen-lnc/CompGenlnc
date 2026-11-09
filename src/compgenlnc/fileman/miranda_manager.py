from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from compgenlnc.consts.fields import PAIR_FIELDS
from compgenlnc.consts.params import LNCRNA, MIRNA
from compgenlnc.consts.paths import BINDING_PREFIX
from compgenlnc.consts.regex import MIRANDA_INFO
from compgenlnc.structs.binding_zone import BindingZone
from compgenlnc.structs.seq_ss_record import SeqSSRecord
from compgenlnc.typing.numpy_dtypes import InteractionPair

if TYPE_CHECKING:
    from _typeshed import StrPath


def read_binding_zone(
    binding_dir: StrPath, interaction: InteractionPair
) -> BindingZone:
    binding_dir = Path(binding_dir).resolve()
    lnc, mir = interaction[PAIR_FIELDS]
    filename = binding_dir / f"{BINDING_PREFIX}{lnc}_{mir}.dat"
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
    binding_dir: StrPath, interactions: Iterable[InteractionPair]
) -> Iterator[BindingZone]:
    return (
        read_binding_zone(binding_dir, interaction)
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
    records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    binding_zones: Sequence[BindingZone],
    *,
    rna: str = LNCRNA,
) -> Iterator[SeqSSRecord]:
    return (
        get_subrecord(*args, rna=rna)
        for args in zip(records, interactions, binding_zones)
    )


def get_subrecord_list_from_dir(
    records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    binding_dir: StrPath,
    *,
    rna: str = LNCRNA,
) -> Iterator[SeqSSRecord]:
    binding_zones = iter_binding_zones(binding_dir, interactions)
    return get_subrecord_list(records, interactions, binding_zones, rna=rna)
