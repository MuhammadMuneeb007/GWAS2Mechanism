"""Thread accounting: SLURM awareness and BLAS oversubscription control."""

from __future__ import annotations

import contextlib
import logging
import os
from collections.abc import Iterator

log = logging.getLogger(__name__)

_BLAS_ENV = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def available_cpus() -> int:
    """CPUs this process may use: SLURM allocation first, then affinity, then count."""
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit() and int(slurm) > 0:
        return int(slurm)
    try:
        return len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def resolve_threads(value: int | str | None) -> int:
    if value in (None, "auto"):
        return max(1, available_cpus())
    return max(1, int(value))


def thread_env(n: int) -> dict[str, str]:
    """Environment for a child process allowed to use ``n`` threads."""
    env = {k: str(n) for k in _BLAS_ENV}
    env["POLARS_MAX_THREADS"] = str(n)
    return env


@contextlib.contextmanager
def limit_blas(n: int) -> Iterator[None]:
    """Cap BLAS/OpenMP threads in-process (prevents oversubscription in pools)."""
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:  # pragma: no cover
        yield
        return
    with threadpool_limits(limits=n):
        yield


def split_jobs(total_threads: int, n_tasks: int, per_job: int = 1) -> tuple[int, int]:
    """Return (n_parallel_jobs, threads_per_job) for a pool of ``n_tasks``."""
    per_job = max(1, per_job)
    jobs = max(1, min(n_tasks, total_threads // per_job))
    return jobs, max(1, total_threads // jobs)
