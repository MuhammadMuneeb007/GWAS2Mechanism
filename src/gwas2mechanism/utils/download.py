"""Verified, resumable and race-safe downloads with bounded retry/fallback."""

from __future__ import annotations

import base64
import hashlib
import logging
import os
import shutil
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import httpx
from tenacity import Retrying, retry_if_exception_type, stop_after_attempt, wait_exponential

from gwas2mechanism.utils.io import atomic_write_json, read_json
from gwas2mechanism.utils.locking import resource_lock
from gwas2mechanism.utils.proc import ToolError, run
from gwas2mechanism.utils.progress import progress

log = logging.getLogger(__name__)

CHUNK = 8 * 1024 * 1024
USER_AGENT = "GWAS2Mechanism/0.1 (+https://github.com/MuhammadMuneeb007/GWAS2Mechanism)"


@dataclass(frozen=True)
class DownloadOptions:
    connections_per_file: int = 4
    min_split_size: str = "16M"
    retry_attempts: int = 8
    retry_initial_seconds: float = 2
    retry_max_seconds: float = 60
    aria2_max_tries: int = 5
    aria2_retry_wait_seconds: int = 5
    timeout_seconds: float = 300
    connect_timeout_seconds: float = 30

    @classmethod
    def from_config(cls, value: object | None) -> DownloadOptions:
        if value is None:
            return cls()
        return cls(**{name: getattr(value, name) for name in cls.__dataclass_fields__})


@dataclass
class DownloadResult:
    url: str
    path: Path
    size: int
    md5: str | None
    status: str
    resumed: bool = False


class DownloadStats:
    """Thread-safe aggregate used by setup performance reporting."""

    def __init__(self) -> None:
        self.files_downloaded = 0
        self.bytes_downloaded = 0
        self.cache_hits = 0
        self.resumed = 0
        self._lock = threading.Lock()

    def observe(self, result: DownloadResult) -> None:
        with self._lock:
            if result.status == "EXISTS":
                self.cache_hits += 1
            else:
                self.files_downloaded += 1
                self.bytes_downloaded += result.size
            if result.resumed:
                self.resumed += 1

    def as_dict(self) -> dict[str, int]:
        return {
            "files_downloaded": self.files_downloaded,
            "bytes_downloaded": self.bytes_downloaded,
            "cache_hits": self.cache_hits,
            "downloads_resumed": self.resumed,
            "downloads_newly_completed": self.files_downloaded,
        }


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    size = path.stat().st_size
    bar = progress(total=size, unit="B", unit_scale=True, unit_divisor=1024, desc=f"Verify {path.name}", disable=size < 512 * 1024 * 1024)
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(CHUNK), b""):
                digest.update(block)
                bar.update(len(block))
    finally:
        bar.close()
    return digest.hexdigest()


def normalize_md5(value: str | None) -> str | None:
    if not value:
        return None
    value = value.strip()
    if len(value) == 32 and all(c in "0123456789abcdefABCDEF" for c in value):
        return value.lower()
    try:
        return base64.b64decode(value).hex()
    except (ValueError, TypeError):
        return None


def _stamp(path: Path) -> Path:
    return path.with_name(path.name + ".verified.json")


def is_verified(path: Path, expected_size: int | None = None, md5: str | None = None) -> bool:
    if not path.is_file():
        return False
    stat = path.stat()
    if expected_size is not None and stat.st_size != expected_size:
        return False
    stamp = _stamp(path)
    if stamp.exists():
        info = read_json(stamp)
        same_file = info.get("size") == stat.st_size and abs(info.get("mtime", 0) - stat.st_mtime) < 1
        if same_file and (md5 is None or info.get("md5") == normalize_md5(md5)):
            return True
    if md5 is not None:
        return md5_file(path) == normalize_md5(md5)
    return expected_size is not None


def write_stamp(path: Path, url: str, md5: str | None = None) -> None:
    stat = path.stat()
    atomic_write_json(
        _stamp(path),
        {"url": url, "size": stat.st_size, "mtime": stat.st_mtime, "md5": normalize_md5(md5) if md5 else None, "verified_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
    )


def is_local(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.scheme in ("", "file") or (len(parsed.scheme) == 1 and os.name == "nt")


def local_path(url: str) -> Path:
    parsed = urlparse(url)
    if parsed.scheme == "file":
        return Path(parsed.path.lstrip("/") if os.name == "nt" else parsed.path)
    return Path(url)


def remote_size(client: httpx.Client, url: str) -> int | None:
    try:
        head = client.head(url)
        if head.status_code == 200 and head.headers.get("content-length"):
            return int(head.headers["content-length"])
        with client.stream("GET", url, headers={"Range": "bytes=0-0"}) as response:
            total = response.headers.get("content-range", "").rsplit("/", 1)[-1]
            return int(total) if total.isdigit() else None
    except httpx.HTTPError:
        return None


class Downloader:
    def __init__(self, use_aria2: list[str] | None = None, timeout: float | None = None, *, options: DownloadOptions | object | None = None, stats: DownloadStats | None = None):
        self.aria2 = use_aria2
        self.options = options if isinstance(options, DownloadOptions) else DownloadOptions.from_config(options)
        if timeout is not None:
            self.options = DownloadOptions(**{**self.options.__dict__, "timeout_seconds": timeout})
        self.stats = stats

    def _client(self) -> httpx.Client:
        return httpx.Client(timeout=httpx.Timeout(self.options.timeout_seconds, connect=self.options.connect_timeout_seconds), headers={"User-Agent": USER_AGENT}, follow_redirects=True)

    def fetch(self, url: str, dest: Path, *, md5: str | None = None, expected_size: int | None = None) -> DownloadResult:
        dest.parent.mkdir(parents=True, exist_ok=True)
        with resource_lock(dest):
            result = self._fetch_locked(url, dest, md5=md5, expected_size=expected_size)
        if self.stats:
            self.stats.observe(result)
        return result

    def _fetch_locked(self, url: str, dest: Path, *, md5: str | None, expected_size: int | None) -> DownloadResult:
        if is_local(url):
            src = local_path(url)
            if not src.exists():
                raise FileNotFoundError(src)
            if is_verified(dest, src.stat().st_size, md5):
                return DownloadResult(url, dest, dest.stat().st_size, md5, "EXISTS")
            part = dest.with_name(dest.name + ".part")
            shutil.copyfile(src, part)
            os.replace(part, dest)
            write_stamp(dest, url, md5)
            return DownloadResult(url, dest, dest.stat().st_size, md5, "COPIED")

        with self._client() as client:
            expected_size = expected_size if expected_size is not None else remote_size(client, url)
            if is_verified(dest, expected_size, md5):
                if not _stamp(dest).exists():
                    write_stamp(dest, url, md5)
                log.info("[cache] %s", dest)
                return DownloadResult(url, dest, dest.stat().st_size, md5, "EXISTS")
            if dest.exists():
                log.warning("Existing %s failed verification; re-downloading", dest)
                dest.unlink()
            part = dest.with_name(dest.name + ".part")
            resumed = part.exists() and part.stat().st_size > 0
            if self.aria2:
                try:
                    self._aria2_with_retry(url, part)
                except (ToolError, OSError) as exc:
                    log.warning("aria2 failed repeatedly for %s; resuming with HTTPX: %s", url, exc)
                    self._httpx_with_retry(client, url, part, expected_size)
            else:
                self._httpx_with_retry(client, url, part, expected_size)

        if not part.exists():
            raise OSError(f"Download produced no file: {url}")
        size = part.stat().st_size
        if expected_size is not None and size != expected_size:
            raise OSError(f"Incomplete download {url}: {size} != {expected_size} bytes")
        if md5 is not None and md5_file(part) != normalize_md5(md5):
            part.unlink()
            raise OSError(f"MD5 mismatch for {url}; corrupt partial file was removed")
        os.replace(part, dest)
        write_stamp(dest, url, md5)
        return DownloadResult(url, dest, size, md5, "DOWNLOADED", resumed=resumed)

    def _aria2_with_retry(self, url: str, part: Path) -> None:
        attempts = max(1, self.options.retry_attempts)
        for attempt in range(1, attempts + 1):
            try:
                self._aria2(url, part)
                return
            except (ToolError, OSError):
                if attempt == attempts:
                    raise
                delay = min(self.options.retry_max_seconds, self.options.retry_initial_seconds * (2 ** (attempt - 1)))
                log.warning("aria2 attempt %d/%d failed for %s; retrying in %.1fs", attempt, attempts, part.name, delay)
                time.sleep(delay)

    def _aria2(self, url: str, part: Path) -> None:
        assert self.aria2 is not None
        connections = str(self.options.connections_per_file)
        run([
            *self.aria2, "--continue=true", f"--max-connection-per-server={connections}", f"--split={connections}",
            f"--min-split-size={self.options.min_split_size}", "--file-allocation=none", "--auto-file-renaming=false",
            "--allow-overwrite=true", f"--max-tries={self.options.aria2_max_tries}",
            f"--retry-wait={self.options.aria2_retry_wait_seconds}", f"--timeout={int(self.options.timeout_seconds)}",
            f"--connect-timeout={int(self.options.connect_timeout_seconds)}", "--console-log-level=warn", "--summary-interval=5",
            "--dir", str(part.parent), "--out", part.name, url,
        ])

    def _httpx_with_retry(self, client: httpx.Client, url: str, part: Path, expected: int | None) -> None:
        retryer = Retrying(
            retry=retry_if_exception_type((httpx.HTTPError, OSError)),
            stop=stop_after_attempt(self.options.retry_attempts),
            wait=wait_exponential(multiplier=self.options.retry_initial_seconds, max=self.options.retry_max_seconds),
            reraise=True,
        )
        retryer(self._httpx, client, url, part, expected)

    def _httpx(self, client: httpx.Client, url: str, part: Path, expected: int | None) -> None:
        existing = part.stat().st_size if part.exists() else 0
        if expected is not None and existing > expected:
            part.unlink()
            existing = 0
        if expected is not None and existing == expected:
            return
        headers = {"Range": f"bytes={existing}-"} if existing else {}
        with client.stream("GET", url, headers=headers) as response:
            if existing and response.status_code == 200:
                existing = 0
            response.raise_for_status()
            total = expected
            if total is None and response.headers.get("content-length"):
                total = existing + int(response.headers["content-length"])
            bar = progress(total=total, initial=existing, unit="B", unit_scale=True, unit_divisor=1024, desc=part.name)
            try:
                with part.open("ab" if existing else "wb") as handle:
                    for block in response.iter_bytes(CHUNK):
                        handle.write(block)
                        bar.update(len(block))
            finally:
                bar.close()


def fetch_range(url: str, n_bytes: int, timeout: float = 60.0) -> bytes:
    """Bounded prefix read (HTTP Range) used for header probing."""
    if is_local(url):
        with local_path(url).open("rb") as handle:
            return handle.read(n_bytes)
    with httpx.Client(timeout=timeout, headers={"User-Agent": USER_AGENT}, follow_redirects=True) as client:
        response = client.get(url, headers={"Range": f"bytes=0-{n_bytes - 1}"})
        response.raise_for_status()
        return response.content[:n_bytes]
