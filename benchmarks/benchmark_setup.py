#!/usr/bin/env python3
"""Compare legacy serial and bounded-parallel setup scheduling with mock files."""

from __future__ import annotations

import argparse
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def mock_resource(path: Path, delay: float) -> None:
    time.sleep(delay)
    path.write_bytes(b"mock-resource")


def run_serial(root: Path, files: int, delay: float) -> float:
    started = time.perf_counter()
    for index in range(files):
        mock_resource(root / f"serial-{index}", delay)
    return time.perf_counter() - started


def run_parallel(root: Path, files: int, delay: float, workers: int) -> float:
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(lambda index: mock_resource(root / f"parallel-{index}", delay), range(files)))
    return time.perf_counter() - started


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", type=int, default=12)
    parser.add_argument("--delay", type=float, default=0.05, help="Mock latency per file in seconds")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="gwas2m-benchmark-") as temporary:
        root = Path(temporary)
        serial = run_serial(root, args.files, args.delay)
        parallel = run_parallel(root, args.files, args.delay, args.workers)
    print(f"Mock resources: {args.files}")
    print(f"Legacy serial: {serial:.3f}s")
    print(f"Bounded parallel ({args.workers} workers): {parallel:.3f}s")
    print(f"Scheduling speedup: {serial / parallel:.2f}x")


if __name__ == "__main__":
    main()
