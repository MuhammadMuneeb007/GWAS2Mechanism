"""Verified, resumable, cache-once downloads.

* existing files are re-used when a verification stamp matches (size + mtime,
  and checksum when one is known) - nothing is re-hashed or re-downloaded;
* ``aria2c`` is used for large parallel/resumable transfers when available;
* otherwise ``httpx`` streams into ``<file>.part`` with HTTP Range resume;
* completion is an atomic ``.part -> final`` rename after size/MD5 checks.
"""

from __future__ import annotations

import base64
import hashlib
import logging
import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from gwas2mechanism.utils.io import atomic_write_json, read_json
from gwas2mechanism.utils.proc import run

log = logging.getLogger(__name__)

CHUNK = 8 * 1024 * 1024
USER_AGENT = "GWAS2Mechanism/0.1 (+https://github.com/MuhammadMuneeb007/GWAS2Mechanism)"


@dataclass
class DownloadResult:
    url: str
    path: Path
    size: int
    md5: str | None
    status: str  # EXISTS | DOWNLOADED | COPIED


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_md5(value: str | None) -> str | None:
    """Accept hex or base64 (Google Cloud Storage style) MD5 values."""
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
    if not path.exists():
        return False
    stamp = _stamp(path)
    stat = path.stat()
    if expected_size is not None and stat.st_size != expected_size:
        return False
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
        {
            "url": url,
            "size": stat.st_size,
            "mtime": stat.st_mtime,
            "md5": normalize_md5(md5) if md5 else None,
            "verified_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
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
        head = client.head(url, follow_redirects=True)
        if head.status_code == 200 and head.headers.get("content-length"):
            return int(head.headers["content-length"])
        with client.stream("GET", url, headers={"Range": "bytes=0-0"}, follow_redirects=True) as r:
            content_range = r.headers.get("content-range", "")
            total = content_range.rsplit("/", 1)[-1]
            return int(total) if total.isdigit() else None
    except httpx.HTTPError:
        return None


class Downloader:
    def __init__(self, use_aria2: list[str] | None = None, timeout: float = 300.0):
        self.aria2 = use_aria2
        self.timeout = timeout

    def _client(self) -> httpx.Client:
        return httpx.Client(
            timeout=httpx.Timeout(self.timeout, connect=30.0),
            headers={"User-Agent": USER_AGENT},
            follow_redirects=True,
        )

    def fetch(
        self,
        url: str,
        dest: Path,
        *,
        md5: str | None = None,
        expected_size: int | None = None,
    ) -> DownloadResult:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if is_local(url):
            src = local_path(url)
            if not src.exists():
                raise FileNotFoundError(src)
            if not is_verified(dest, src.stat().st_size, md5):
                tmp = dest.with_name(dest.name + ".part")
                shutil.copyfile(src, tmp)
                os.replace(tmp, dest)
                write_stamp(dest, url, md5)
                return DownloadResult(url, dest, dest.stat().st_size, md5, "COPIED")
            return DownloadResult(url, dest, dest.stat().st_size, md5, "EXISTS")

        with self._client() as client:
            if expected_size is None:
                expected_size = remote_size(client, url)
            if is_verified(dest, expected_size, md5):
                if not _stamp(dest).exists():
                    write_stamp(dest, url, md5)
                log.info("[cache] %s", dest)
                return DownloadResult(url, dest, dest.stat().st_size, md5, "EXISTS")
            if dest.exists():
                log.warning("Existing %s failed verification; re-downloading", dest)
                dest.unlink()
            part = dest.with_name(dest.name + ".part")
            if self.aria2:
                self._aria2(url, part)
            else:
                self._httpx(client, url, part, expected_size)
        size = part.stat().st_size
        if expected_size is not None and size != expected_size:
            raise OSError(f"Incomplete download {url}: {size} != {expected_size} bytes")
        if md5 is not None and md5_file(part) != normalize_md5(md5):
            part.unlink()
            raise OSError(f"MD5 mismatch for {url}")
        os.replace(part, dest)
        write_stamp(dest, url, md5)
        return DownloadResult(url, dest, size, md5, "DOWNLOADED")

    def _aria2(self, url: str, part: Path) -> None:
        assert self.aria2 is not None
        run(
            [
                *self.aria2,
                "--continue=true",
                "--max-connection-per-server=8",
                "--split=8",
                "--min-split-size=16M",
                "--auto-file-renaming=false",
                "--allow-overwrite=true",
                "--console-log-level=warn",
                "--dir",
                str(part.parent),
                "--out",
                part.name,
                url,
            ]
        )

    @retry(
        retry=retry_if_exception_type((httpx.HTTPError, OSError)),
        stop=stop_after_attempt(8),
        wait=wait_exponential(multiplier=2, max=60),
        reraise=True,
    )
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
                existing = 0  # server ignored Range; restart
            response.raise_for_status()
            with open(part, "ab" if existing else "wb") as handle:
                for block in response.iter_bytes(CHUNK):
                    handle.write(block)


def fetch_range(url: str, n_bytes: int, timeout: float = 60.0) -> bytes:
    """Bounded prefix read (HTTP Range) used for header probing."""
    if is_local(url):
        with open(local_path(url), "rb") as handle:
            return handle.read(n_bytes)
    with httpx.Client(
        timeout=timeout, headers={"User-Agent": USER_AGENT}, follow_redirects=True
    ) as client:
        response = client.get(url, headers={"Range": f"bytes=0-{n_bytes - 1}"})
        response.raise_for_status()
        return response.content[:n_bytes]
