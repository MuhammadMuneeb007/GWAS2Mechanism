"""Population classification from GWAS Catalog ancestry metadata.

Key rules
---------
* A study's population is derived from its **discovery (initial) stage**.
* One 1000 Genomes super-population -> that label (``approximate`` flagged when
  the Catalog category only loosely corresponds, e.g. admixed groups).
* Two or more super-populations -> ``MULTI``; such a study is analysed as
  ``COMBINED_PUBLISHED`` and is **never split** into per-ancestry files.
* Categories without a 1000 Genomes counterpart -> ``UNMAPPED``; missing or
  "NR" metadata -> ``UNKNOWN``. Neither is forced into a super-population.
* The original Catalog category text is always retained.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

from gwas2mechanism.config import PopulationMapEntry
from gwas2mechanism.constants import (
    COMBINED_PUBLISHED,
    MULTI,
    SUPERPOPULATIONS,
    UNKNOWN_POPULATION,
    UNMAPPED_POPULATION,
)
from gwas2mechanism.discovery.models import AncestryRecord, CandidateStudy

#: GWAS Catalog broad ancestral categories (Morales et al. 2018 framework).
KNOWN_CATEGORIES: tuple[str, ...] = (
    "Greater Middle Eastern (Middle Eastern, North African or Persian)",
    "African American or Afro-Caribbean",
    "Hispanic or Latin American",
    "Other admixed ancestry",
    "Aboriginal Australian",
    "African unspecified",
    "Sub-Saharan African",
    "South East Asian",
    "Native American",
    "East Asian",
    "South Asian",
    "Asian unspecified",
    "European",
    "Oceanian",
    "Other",
    "NR",
)

_UNKNOWN_CATEGORIES = {"nr", "not reported", "other", "unknown", ""}
_COUNT_PREFIX = re.compile(r"^\s*(\d[\d,]*)\s+(.*)$")


def parse_ancestry_string(text: str, stage: str = "initial") -> list[AncestryRecord]:
    """Parse REST v2 ``discovery_ancestry`` entries like ``437667 European (U.S.)``.

    Multiple categories may be comma-joined for one sample block
    (e.g. ``1200 European, East Asian``); each is kept, with the count only on
    the block (not duplicated).
    """
    text = (text or "").strip()
    n: int | None = None
    match = _COUNT_PREFIX.match(text)
    if match:
        n = int(match.group(1).replace(",", ""))
        text = match.group(2)
    records: list[AncestryRecord] = []
    remaining = text
    while remaining:
        remaining = remaining.lstrip(" ,;")
        if not remaining:
            break
        category = next((c for c in KNOWN_CATEGORIES if remaining.startswith(c)), None)
        if category is None:
            head, _, tail = remaining.partition(",")
            head = head.strip()
            country = None
            paren = re.match(r"^(.*?)\s*\((.*)\)\s*$", head)
            if paren:
                head, country = paren.group(1).strip(), paren.group(2)
            records.append(AncestryRecord(stage=stage, category=head, country=country))
            remaining = tail
            continue
        remaining = remaining[len(category) :].strip()
        country = None
        if remaining.startswith("("):
            close = remaining.find(")")
            country = remaining[1:close] if close > 0 else remaining[1:]
            remaining = remaining[close + 1 :] if close > 0 else ""
        records.append(AncestryRecord(stage=stage, category=category, country=country))
    if not records:
        records = [AncestryRecord(stage=stage, category="NR")]
    if n is not None:
        records[0].n = n
    return records


def classify_category(category: str, mapping: dict[str, PopulationMapEntry]) -> tuple[str, bool]:
    key = category.strip()
    if key.lower() in _UNKNOWN_CATEGORIES:
        return UNKNOWN_POPULATION, False
    entry = mapping.get(key)
    if entry is None:
        lowered = {k.lower(): v for k, v in mapping.items()}
        entry = lowered.get(key.lower())
    if entry is None:
        return UNMAPPED_POPULATION, False
    if entry.population not in SUPERPOPULATIONS:
        raise ValueError(f"population mapping for {category!r} targets unknown label {entry.population!r}")
    return entry.population, entry.approximate


def classify_records(records: Iterable[AncestryRecord], mapping: dict[str, PopulationMapEntry]) -> list[AncestryRecord]:
    out = []
    for record in records:
        record.population, record.approximate = classify_category(record.category, mapping)
        out.append(record)
    return out


def classify_study(study: CandidateStudy, mapping: dict[str, PopulationMapEntry]) -> CandidateStudy:
    records = classify_records(study.ancestries, mapping)
    discovery = [r for r in records if r.stage == "initial"] or records
    study.reported_ancestry = "; ".join(
        dict.fromkeys(f"{r.category}{f' ({r.country})' if r.country else ''}" for r in discovery)
    ) or None

    labels = {r.population for r in discovery}
    known = labels - {UNKNOWN_POPULATION}
    if not discovery or not known:
        study.population = UNKNOWN_POPULATION
        study.population_reason = "no usable discovery-stage ancestry metadata"
    elif len(known) == 1:
        (label,) = known
        study.population = label
        study.population_approximate = any(r.approximate for r in discovery if r.population == label)
        if label == UNMAPPED_POPULATION:
            study.population_reason = "ancestry category has no 1000 Genomes super-population counterpart"
        else:
            study.population_reason = (
                f"all discovery samples map to {label}"
                + (" (approximate mapping)" if study.population_approximate else "")
                + (" ; some samples lack ancestry labels" if UNKNOWN_POPULATION in labels else "")
            )
    else:
        study.population = MULTI
        study.population_reason = "discovery sample spans " + ", ".join(sorted(known))

    if study.population in SUPERPOPULATIONS:
        study.analysis_label = study.population
    elif study.population == MULTI:
        # Aggregated multi-ancestry statistics: analysed as-is, never split.
        study.analysis_label = COMBINED_PUBLISHED
    else:
        study.analysis_label = None

    structured = [r.n for r in discovery if r.n is not None]
    if structured and study.n_total is None:
        study.n_total = int(sum(structured))
        study.n_source = "structured_ancestry"
    return study


# ----------------------------------------------------------------------------
# Sample size parsing - never fabricated: unparseable -> None
# ----------------------------------------------------------------------------

_CASES = re.compile(r"(\d[\d,]*)\s+(?:[^,;]*?\s)?cases?\b", re.IGNORECASE)
_CONTROLS = re.compile(r"(\d[\d,]*)\s+(?:[^,;]*?\s)?controls?\b", re.IGNORECASE)
_INDIVIDUALS = re.compile(
    r"(\d[\d,]*)\s+(?:[^,;]*?\s)?(?:individuals?|participants?|people|subjects?|samples?|women|men|adults|children)\b",
    re.IGNORECASE,
)


def _sum(pattern: re.Pattern[str], text: str) -> int | None:
    values = [int(m.replace(",", "")) for m in pattern.findall(text)]
    return sum(values) if values else None


def parse_sample_size(text: str | None) -> dict[str, int | None]:
    """Extract cases/controls/total from Catalog ``initial_sample_size`` text."""
    if not text:
        return {"n_cases": None, "n_controls": None, "n_total": None}
    cases = _sum(_CASES, text)
    controls = _sum(_CONTROLS, text)
    individuals = _sum(_INDIVIDUALS, text)
    if cases is not None and controls is not None:
        total = cases + controls + (individuals or 0)
    else:
        total = individuals
    return {"n_cases": cases, "n_controls": controls, "n_total": total}


def apply_sample_size(study: CandidateStudy) -> CandidateStudy:
    parsed = parse_sample_size(study.initial_sample_size_text)
    study.n_cases = parsed["n_cases"]
    study.n_controls = parsed["n_controls"]
    if study.n_total is None and parsed["n_total"] is not None:
        study.n_total = parsed["n_total"]
        study.n_source = "text_parsed"
    if study.n_cases is not None and study.n_controls is not None:
        study.trait_type = "binary"
    elif study.n_total is not None and study.n_cases is None:
        study.trait_type = "quantitative_or_unspecified"
    return study
