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
    """Read the binding zone of an interaction.

    Parameters
    ----------
    binding_dir : str or os.PathLike
        The directory containing the binding zone file.
    interaction : InteractionPair
        The interaction for which to read the binding zone.

    Returns
    -------
    BindingZone
        An object containing the ranges where the lncRNA and miRNA
        interact.
    """
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


def read_binding_zones_from_list(
    binding_dir: StrPath, interactions: Iterable[InteractionPair]
) -> Iterator[BindingZone]:
    """Read the binding zones from an interaction set.

    Parameters
    ----------
    binding_dir : str or os.PathLike
        The directory containing the binding zone files.
    interactions : iterable object of InteractionPair
        Group of interactions to read the binding zones from.

    Returns
    -------
    Generator of BindingZone
        An iterator of BindingZone containing the ranges where the
        lncRNA and miRNA interact in each interaction.
    """
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
    """Get the subrecord for an RNA in an interaction.

    Parameters
    ----------
    record : SeqSSRecord
        The RNA record to get the subrecord from.
    interaction : InteractionPair
        The interaction pair to identify the subrecord.
    binding_zone : BindingZone
        The binding zone to identify the subrecord.
    rna : str, optional
        The RNA type to identify the subrecord. Default: `LNCRNA`.

    Returns
    -------
    SeqSSRecord
        The subrecord for the RNA in the interaction.
    """
    if rna == LNCRNA:
        start = binding_zone.lnc_start - 1
        stop = binding_zone.lnc_end
    elif rna == MIRNA:
        start = binding_zone.mir_start - 1
        stop = binding_zone.mir_end
    return SeqSSRecord(
        "_".join(interaction[PAIR_FIELDS]),
        record.seq[start:stop],
        record.ss[start:stop],
    )


def get_subrecord_list(
    records: Sequence[SeqSSRecord],
    interactions: Sequence[InteractionPair],
    binding_zones: Sequence[BindingZone],
    *,
    rna: str = LNCRNA,
) -> Iterator[SeqSSRecord]:
    """Get the subrecord for an RNA in an interaction.

    Parameters
    ----------
    records : sequence of SeqSSRecord
        Sequence of RNA records to get the subrecords from.
    interactions : sequence of InteractionPair
        Sequence of interaction pairs to identify the subrecords.
    binding_zones : sequence of BindingZone
        Sequence of binding zones to identify the subrecords.
    rna : str, optional
        The RNA type to identify the subrecords. Default: `LNCRNA`.

    Returns
    -------
    Generator of SeqSSRecord
        An iterator that yields the subrecords for the RNAs in the
        interactions.
    """
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
    """Get the subrecord for an RNA in an interaction.

    Parameters
    ----------
    records : sequence of SeqSSRecord
        Sequence of RNA records to get the subrecords from.
    interactions : sequence of InteractionPair
        Sequence of interaction pairs to identify the subrecords.
    binding_dir : str or os.PathLike
        Path of the directory containing the binding zone files.
    rna : str, optional
        The RNA type to identify the subrecords. Default: `LNCRNA`.

    Returns
    -------
    Generator of SeqSSRecord
        An iterator that yields the subrecords for the RNAs in the
        interactions.
    """
    binding_zones = read_binding_zones_from_list(binding_dir, interactions)
    return get_subrecord_list(records, interactions, binding_zones, rna=rna)
