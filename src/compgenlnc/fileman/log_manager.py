from collections.abc import Iterable, Iterator
from typing import TypeVar

from tqdm import tqdm

from compgenlnc.config.var_config import LOG_FILE

T = TypeVar("T")


def send_msg(
    msg: str, flush: bool = True, seek: int | None = None, seek_mode: int = 0
) -> None:
    if LOG_FILE is None:
        return
    if seek is not None:
        LOG_FILE.seek(seek, seek_mode)
    LOG_FILE.write(msg)
    if flush:
        LOG_FILE.flush()


def track_progress(
    it: Iterable[T] | None = None,
    desc: str | None = None,
    total: int | None = None,
    unit: str = "it",
    initial: int = 0,
    steps: int = 1,
) -> Iterator[T]:
    if LOG_FILE is None:
        return it
    progress = tqdm(
        desc=desc, total=total, file=LOG_FILE, unit=unit, initial=initial
    )

    def iterator():
        for elem in it:
            yield elem
            progress.update(min(steps, total - progress.n))

    return (x for x in iterator())
