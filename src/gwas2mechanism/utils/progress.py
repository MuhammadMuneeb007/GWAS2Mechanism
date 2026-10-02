"""TTY- and SLURM-friendly tqdm helpers."""

from __future__ import annotations

import os
import sys
from collections.abc import Iterable, Iterator
from typing import TypeVar

from tqdm.auto import tqdm

T = TypeVar("T")


def is_batch() -> bool:
    return bool(os.environ.get("SLURM_JOB_ID")) or not sys.stderr.isatty()


def progress(iterable: Iterable[T] | None = None, **kwargs) -> tqdm | Iterator[T]:
    """Create a progress bar without flooding redirected/SLURM logs."""
    kwargs.setdefault("dynamic_ncols", not is_batch())
    kwargs.setdefault("mininterval", 5.0 if is_batch() else 0.1)
    kwargs.setdefault("maxinterval", 30.0)
    kwargs.setdefault("leave", True)
    return tqdm(iterable, **kwargs)
