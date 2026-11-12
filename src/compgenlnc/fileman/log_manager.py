from collections.abc import Iterable, Iterator
from typing import TypeVar

from tqdm import tqdm

import compgenlnc.config as config

T = TypeVar("T")


def send_msg(
    msg: str, /, flush: bool = True, seek: int = 0, seek_mode: int = 0
) -> None:
    """Write a message to the log file.

    Parameters
    ----------
    msg : str
        The message to write to the log file.
    flush : bool, optional
        Whether to flush the log file after writing the message.
        Default: `True`.
    seek : int or None, optional
        The position to seek in the log file before writing the message.
        Default: `0`.
    seek_mode : int, optional
        The mode to use when seeking in the log file. Default: `0`.
    """
    if config.log_file is None:
        return
    if seek is not None:
        config.log_file.seek(seek, seek_mode)
    config.log_file.write(msg)
    if flush:
        config.log_file.flush()


def track_progress(
    it: Iterable[T],
    /,
    desc: str | None = None,
    total: int | None = None,
    unit: str = "it",
    initial: int = 0,
    steps: int = 1,
) -> Iterator[T]:
    """Display a progress bar for an iterable in the log file.

    Parameters
    ----------
    it : iterable
        The iterable to track progress for.
    desc : str, optional
        The description of the progress bar. Default: `None`.
    total : int, optional
        The total number of items in the iterable. Default: `None`.
    unit : str, optional
        The unit of the elements in the iterable. Default: `"it"`.
    initial : int, optional
        The initial value of the progress bar. Default: `0`.
    steps : int, optional
        The number of items to update the progress bar by. Default: `1`.

    Returns
    -------
    Generator
        An iterator that yields the elements of the iterable with a
        progress bar.
    """
    if config.log_file is None:
        return it
    if total is None:
        it = list(it)
        total = len(it)
    progress = tqdm(
        desc=desc, total=total, file=config.log_file, unit=unit, initial=initial
    )

    def iterator():
        for elem in it:
            yield elem
            progress.update(min(steps, total - progress.n))

    return (x for x in iterator())
