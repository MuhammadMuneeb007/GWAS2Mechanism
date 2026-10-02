"""Small dependency-free inter-process lock based on atomic directory creation."""

from __future__ import annotations

import json
import os
import socket
import threading
import time
from pathlib import Path

_THREAD_LOCKS: dict[str, threading.RLock] = {}
_THREAD_LOCKS_GUARD = threading.Lock()


def _thread_lock(path: Path) -> threading.RLock:
    key = str(path.resolve())
    with _THREAD_LOCKS_GUARD:
        return _THREAD_LOCKS.setdefault(key, threading.RLock())


class FileLock:
    """Wait for exclusive ownership of ``path`` and clean stale local locks."""

    def __init__(self, path: Path, *, timeout: float = 3600, poll: float = 0.25, stale_after: float = 86400):
        self.path = path
        self.timeout = timeout
        self.poll = poll
        self.stale_after = stale_after
        self._owned = False
        self._local = _thread_lock(path)

    def acquire(self) -> None:
        if not self._local.acquire(timeout=self.timeout):
            raise TimeoutError(f"Timed out waiting for in-process lock {self.path}")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        started = time.monotonic()
        try:
            while True:
                try:
                    self.path.mkdir()
                    self._owned = True
                    (self.path / "owner.json").write_text(
                        json.dumps({"pid": os.getpid(), "hostname": socket.gethostname(), "created": time.time()}),
                        encoding="utf-8",
                    )
                    return
                except FileExistsError:
                    if self._is_stale():
                        self._break_stale()
                        continue
                    if time.monotonic() - started >= self.timeout:
                        raise TimeoutError(f"Timed out waiting for lock {self.path}") from None
                    time.sleep(self.poll)
        except Exception:
            self._local.release()
            raise

    def _is_stale(self) -> bool:
        try:
            info = json.loads((self.path / "owner.json").read_text(encoding="utf-8"))
            age = time.time() - float(info.get("created", 0))
            if info.get("hostname") == socket.gethostname():
                try:
                    os.kill(int(info["pid"]), 0)
                    return False
                except (OSError, KeyError, ValueError):
                    return True
            return age > self.stale_after
        except (OSError, ValueError, json.JSONDecodeError):
            try:
                return time.time() - self.path.stat().st_mtime > self.stale_after
            except OSError:
                return False

    def _break_stale(self) -> None:
        try:
            (self.path / "owner.json").unlink(missing_ok=True)
            self.path.rmdir()
        except OSError:
            pass

    def release(self) -> None:
        if not self._owned:
            return
        try:
            for _ in range(20):
                try:
                    (self.path / "owner.json").unlink(missing_ok=True)
                    self.path.rmdir()
                    break
                except PermissionError:
                    time.sleep(0.05)
        finally:
            self._owned = False
            self._local.release()

    def __enter__(self) -> FileLock:
        self.acquire()
        return self

    def __exit__(self, *_exc) -> None:
        self.release()


def resource_lock(path: Path, **kwargs) -> FileLock:
    return FileLock(path.with_name(path.name + ".lock"), **kwargs)
