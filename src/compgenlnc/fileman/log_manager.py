from collections.abc import Iterable

from tqdm import tqdm

from compgenlnc.config.var_config import LOG_FILE


def send_msg(
    msg: str, flush: bool = True, seek: int | None = None, seek_mode: int = 0
) -> None:
    if seek is not None:
        LOG_FILE.seek(seek, seek_mode)
    LOG_FILE.write(msg)
    if flush:
        LOG_FILE.flush()


def track_progress(
    it: Iterable,
    desc: str | None = None,
    total: int | None = None,
    initial: int = 0,
) -> tqdm:
    return tqdm(it, desc=desc, total=total, file=LOG_FILE, initial=initial)
