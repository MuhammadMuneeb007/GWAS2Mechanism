"""Random access to an (uncompressed) indexed FASTA.

Uses ``pysam.FastaFile`` when importable (handles bgzip), otherwise a small
pure-Python reader driven by the standard ``.fai`` index (built here when
absent). Used for build validation and REF-allele checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class _FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


def build_fai(fasta: Path) -> Path:
    fai = fasta.with_name(fasta.name + ".fai")
    rows: list[str] = []
    with open(fasta, "rb") as handle:
        name = None
        length = offset = line_bases = line_width = 0
        pos = 0
        for raw in handle:
            if raw.startswith(b">"):
                if name is not None:
                    rows.append(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}")
                name = raw[1:].split()[0].decode()
                length = line_bases = line_width = 0
                offset = pos + len(raw)
            else:
                stripped = raw.rstrip(b"\r\n")
                if line_bases == 0:
                    line_bases = len(stripped)
                    line_width = len(raw)
                length += len(stripped)
            pos += len(raw)
        if name is not None:
            rows.append(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_width}")
    fai.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return fai


class FastaReader:
    def __init__(self, path: Path):
        self.path = Path(path)
        self._pysam = None
        try:  # pragma: no cover - exercised only where pysam exists
            import pysam

            self._pysam = pysam.FastaFile(str(self.path))
            self._names = set(self._pysam.references)
            return
        except Exception:
            self._pysam = None
        fai = self.path.with_name(self.path.name + ".fai")
        if not fai.exists():
            build_fai(self.path)
        self._index: dict[str, _FaiEntry] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 5:
                self._index[parts[0]] = _FaiEntry(*(int(x) for x in parts[1:5]))
        self._names = set(self._index)
        self._handle = open(self.path, "rb")  # noqa: SIM115

    def contig(self, chrom: str) -> str | None:
        for candidate in (chrom, f"chr{chrom}", chrom.removeprefix("chr")):
            if candidate in self._names:
                return candidate
        if chrom in ("MT", "M"):
            for candidate in ("chrM", "MT", "M"):
                if candidate in self._names:
                    return candidate
        return None

    def fetch(self, chrom: str, start1: int, length: int) -> str | None:
        """Return ``length`` bases starting at 1-based ``start1``."""
        name = self.contig(str(chrom))
        if name is None or start1 < 1:
            return None
        if self._pysam is not None:  # pragma: no cover
            return self._pysam.fetch(name, start1 - 1, start1 - 1 + length).upper()
        entry = self._index[name]
        if start1 - 1 + length > entry.length:
            return None
        out = []
        pos0 = start1 - 1
        remaining = length
        while remaining > 0:
            line_no, col = divmod(pos0, entry.line_bases)
            take = min(remaining, entry.line_bases - col)
            self._handle.seek(entry.offset + line_no * entry.line_width + col)
            out.append(self._handle.read(take).decode())
            pos0 += take
            remaining -= take
        return "".join(out).upper()

    def close(self) -> None:
        if self._pysam is None and hasattr(self, "_handle"):
            self._handle.close()

    def __enter__(self) -> FastaReader:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
