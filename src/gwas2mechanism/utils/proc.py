"""Subprocess execution with logging and informative failures."""

from __future__ import annotations

import logging
import os
import shlex
import subprocess
from collections.abc import Sequence
from pathlib import Path

log = logging.getLogger(__name__)


class ToolError(RuntimeError):
    """An external command failed; carries the tail of its output."""


def run(
    cmd: Sequence[str | os.PathLike[str]],
    *,
    log_file: Path | None = None,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    check: bool = True,
    capture_tail: int = 60,
) -> subprocess.CompletedProcess[str]:
    argv = [str(c) for c in cmd]
    log.info("$ %s", " ".join(shlex.quote(a) for a in argv))
    full_env = os.environ.copy()
    if env:
        full_env.update({k: str(v) for k, v in env.items()})
    proc = subprocess.run(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=full_env,
        cwd=str(cwd) if cwd else None,
        errors="replace",
    )
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_file.write_text(
            "$ " + " ".join(shlex.quote(a) for a in argv) + "\n" + (proc.stdout or ""),
            encoding="utf-8",
        )
    if check and proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-capture_tail:])
        raise ToolError(
            f"Command failed with exit status {proc.returncode}:\n"
            f"  {' '.join(argv)}\n--- output (tail) ---\n{tail}"
        )
    return proc
