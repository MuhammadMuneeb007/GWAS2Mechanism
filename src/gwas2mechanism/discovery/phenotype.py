"""Phenotype resolution and study-trait matching.

Resolution order (never a bare substring test):

1. user-pinned ontology ids (``phenotype_resolution.efo_ids``);
2. an ontology term whose **label** equals the input (normalised);
3. an ontology term with an **exact synonym** equal to the input;
4. a controlled fuzzy match (flagged, and not selectable unless ``allow_fuzzy``);
5. otherwise the input is used as free text for exact reported-trait matching.

Related terms returned by ontology search (subtypes, "age of onset of ...",
medication-use measurements, ...) are *recorded as rejected with a reason*
rather than silently included. Descendant terms are only added when
``include_descendants`` is enabled.
"""

from __future__ import annotations

import difflib
import logging
import re
import urllib.parse
from collections.abc import Callable, Iterable
from typing import Any

import httpx

from gwas2mechanism.config import PhenotypeResolutionConfig
from gwas2mechanism.discovery.models import CandidateStudy, OntologyTerm, PhenotypeResolution
from gwas2mechanism.paths import slugify

log = logging.getLogger(__name__)

# match_type -> score. Documented, interpretable, and exported per study.
MATCH_SCORES: dict[str, float] = {
    "exact_mapped_trait": 1.0,
    "exact_reported_trait": 0.95,
    "synonym_mapped_trait": 0.9,
    "synonym_reported_trait": 0.85,
    "descendant_term": 0.7,
    "fuzzy_reported_trait": 0.5,
    "none": 0.0,
}

_QUALIFIER = re.compile(r"\s*[\(\[][^\)\]]*[\)\]]\s*")
_STRATIFIED = re.compile(
    r"\b(males?|females?|men|women|sex[- ]stratified|age[- ]stratified|early[- ]onset|late[- ]onset)\b",
    re.IGNORECASE,
)


def normalize_text(text: str | None) -> str:
    if not text:
        return ""
    text = text.lower().replace("_", " ")
    text = re.sub(r"[^\w\s'-]", " ", text)
    return " ".join(text.split())


def normalize_reported(text: str | None) -> tuple[str, str | None]:
    """Strip parenthetical qualifiers: ``Migraine (PheCode 340)`` -> ``migraine``."""
    if not text:
        return "", None
    qualifiers = " ".join(q.strip(" ()[]") for q in _QUALIFIER.findall(text)) or None
    return normalize_text(_QUALIFIER.sub(" ", text)), qualifiers


def is_stratified_subset(reported: str | None) -> bool:
    return bool(reported and _STRATIFIED.search(reported))


def canonical_term_id(value: str) -> str:
    """``MONDO:0005277`` / IRI / ``MONDO_0005277`` -> ``MONDO_0005277``."""
    value = value.strip().rsplit("/", 1)[-1]
    return value.replace(":", "_")


# ----------------------------------------------------------------------------
# Ontology lookup (OLS4) - injectable for tests / offline use
# ----------------------------------------------------------------------------

OntologySearch = Callable[[str, bool], list[OntologyTerm]]


def ols_search_factory(base_url: str, timeout: float = 30.0) -> OntologySearch:
    def search(query: str, exact: bool) -> list[OntologyTerm]:
        params = {
            "q": query,
            "ontology": "efo",
            "type": "class",
            "rows": "25",
            "queryFields": "label,synonym",
            "fieldList": "iri,label,short_form,obo_id,synonym",
        }
        if exact:
            params["exact"] = "true"
        try:
            response = httpx.get(f"{base_url}/search", params=params, timeout=timeout)
            response.raise_for_status()
            docs = response.json().get("response", {}).get("docs", [])
        except (httpx.HTTPError, ValueError) as exc:
            log.warning("OLS search failed for %r: %s", query, exc)
            return []
        terms = []
        for doc in docs:
            if not doc.get("short_form") or not doc.get("label"):
                continue
            terms.append(
                OntologyTerm(
                    term_id=canonical_term_id(doc["short_form"]),
                    label=doc["label"],
                    synonyms=tuple(doc.get("synonym") or ()),
                    iri=doc.get("iri"),
                    source="ols",
                )
            )
        return terms

    return search


def ols_descendants_factory(base_url: str, timeout: float = 30.0) -> Callable[[OntologyTerm], list[str]]:
    def descendants(term: OntologyTerm) -> list[str]:
        if not term.iri:
            return []
        iri = urllib.parse.quote(urllib.parse.quote(term.iri, safe=""), safe="")
        url = f"{base_url}/ontologies/efo/terms/{iri}/hierarchicalDescendants"
        ids: list[str] = []
        try:
            while url:
                response = httpx.get(url, params={"size": 500}, timeout=timeout)
                response.raise_for_status()
                payload = response.json()
                for item in payload.get("_embedded", {}).get("terms", []):
                    if item.get("short_form"):
                        ids.append(canonical_term_id(item["short_form"]))
                url = payload.get("_links", {}).get("next", {}).get("href")
        except (httpx.HTTPError, ValueError) as exc:
            log.warning("OLS descendant lookup failed for %s: %s", term.term_id, exc)
        return ids

    return descendants


# ----------------------------------------------------------------------------
# Resolver
# ----------------------------------------------------------------------------


class PhenotypeResolver:
    def __init__(
        self,
        config: PhenotypeResolutionConfig,
        ontology_search: OntologySearch | None = None,
        descendants: Callable[[OntologyTerm], list[str]] | None = None,
        extra_terms: Iterable[OntologyTerm] = (),
    ):
        self.config = config
        self.search = ontology_search
        self.descendants = descendants
        self.extra_terms = list(extra_terms)  # e.g. Catalog trait list in bulk mode

    def resolve(self, phenotype: str) -> PhenotypeResolution:
        query = normalize_text(phenotype)
        user_syn = [normalize_text(s) for s in self.config.synonyms if s.strip()]
        candidates: list[OntologyTerm] = []
        if self.search is not None:
            candidates += self.search(phenotype, True)
        candidates += self.extra_terms
        candidates = _dedupe(candidates)

        accepted: list[OntologyTerm] = []
        resolution_type = "free_text"

        pinned = {canonical_term_id(x) for x in self.config.efo_ids}
        if pinned:
            by_id = {t.term_id: t for t in candidates}
            for term_id in sorted(pinned):
                accepted.append(by_id.get(term_id, OntologyTerm(term_id, term_id, source="user")))
            resolution_type = "user_pinned"

        if not accepted:
            exact_label = [t for t in candidates if normalize_text(t.label) == query]
            exact_syn = [
                t
                for t in candidates
                if query in {normalize_text(s) for s in t.synonyms} and t not in exact_label
            ]
            if exact_label:
                accepted, resolution_type = exact_label, "exact_label"
            elif exact_syn:
                accepted, resolution_type = exact_syn, "exact_synonym"

        rejected: list[dict[str, Any]] = []
        if not accepted and self.search is not None:
            related = _dedupe(self.search(phenotype, False) + candidates)
            best, ratio = _best_fuzzy(query, related)
            if best is not None and ratio >= self.config.fuzzy_min_ratio:
                accepted, resolution_type = [best], "fuzzy"
            candidates = related

        accepted_ids = {t.term_id for t in accepted}
        for term in candidates:
            if term.term_id not in accepted_ids:
                rejected.append(
                    {
                        "term_id": term.term_id,
                        "label": term.label,
                        "reason": "related ontology term, not an exact label/synonym match; "
                        "pin via phenotype_resolution.efo_ids or enable include_descendants",
                    }
                )

        canonical = accepted[0].label if accepted else phenotype.strip()
        synonyms = {query, normalize_text(canonical), *user_syn}
        for term in accepted:
            synonyms.add(normalize_text(term.label))
            synonyms.update(normalize_text(s) for s in term.synonyms)
        synonyms.discard("")

        descendant_ids: list[str] = []
        if self.config.include_descendants and self.descendants is not None:
            for term in accepted:
                descendant_ids += self.descendants(term)

        return PhenotypeResolution(
            phenotype_input=phenotype,
            phenotype_canonical=canonical,
            phenotype_slug=slugify(phenotype),
            terms=accepted,
            synonyms=sorted(synonyms),
            descendant_ids=sorted(set(descendant_ids) - accepted_ids),
            rejected_terms=rejected[:200],
            resolution_type=resolution_type,
            source="ols" if self.search is not None else "offline",
        )


def _dedupe(terms: list[OntologyTerm]) -> list[OntologyTerm]:
    seen: dict[str, OntologyTerm] = {}
    for term in terms:
        if term.term_id not in seen or (not seen[term.term_id].synonyms and term.synonyms):
            seen[term.term_id] = term
    return list(seen.values())


def _best_fuzzy(query: str, terms: list[OntologyTerm]) -> tuple[OntologyTerm | None, float]:
    best, best_ratio = None, 0.0
    for term in terms:
        for text in (term.label, *term.synonyms):
            ratio = difflib.SequenceMatcher(None, query, normalize_text(text)).ratio()
            if ratio > best_ratio:
                best, best_ratio = term, ratio
    return best, best_ratio


# ----------------------------------------------------------------------------
# Study matching
# ----------------------------------------------------------------------------


def match_study(study: CandidateStudy, resolution: PhenotypeResolution, fuzzy_min_ratio: float = 0.92) -> CandidateStudy:
    """Annotate ``study`` with match_type / match_score / match_reason (WHY it matched)."""
    ids = set(study.mapped_trait_ids)
    labels = {normalize_text(x): x for x in study.mapped_trait_labels}
    canonical = normalize_text(resolution.phenotype_canonical)
    synonyms = set(resolution.synonyms)
    reported, _ = normalize_reported(study.reported_trait)

    hits: list[tuple[str, str]] = []
    shared = ids & resolution.term_ids
    if shared:
        hits.append(
            (
                "exact_mapped_trait",
                f"mapped ontology term {sorted(shared)[0]} is the resolved phenotype term",
            )
        )
    if canonical in labels:
        hits.append(("exact_mapped_trait", f"mapped trait '{labels[canonical]}' equals canonical phenotype"))
    if reported and reported in (canonical, normalize_text(resolution.phenotype_input)):
        hits.append(("exact_reported_trait", f"reported trait '{study.reported_trait}' equals phenotype"))
    syn_mapped = synonyms & set(labels)
    if syn_mapped:
        hits.append(("synonym_mapped_trait", f"mapped trait '{labels[sorted(syn_mapped)[0]]}' is a synonym"))
    if reported and reported in synonyms:
        hits.append(("synonym_reported_trait", f"reported trait '{study.reported_trait}' is a synonym"))
    desc = ids & set(resolution.descendant_ids)
    if desc:
        hits.append(("descendant_term", f"mapped term {sorted(desc)[0]} is a descendant (include_descendants on)"))
    if not hits and reported:
        ratio = max(
            (difflib.SequenceMatcher(None, reported, s).ratio() for s in synonyms), default=0.0
        )
        if ratio >= fuzzy_min_ratio:
            hits.append(("fuzzy_reported_trait", f"reported trait fuzzy ratio {ratio:.2f}"))

    if hits:
        best = max(hits, key=lambda h: MATCH_SCORES[h[0]])
        study.match_type, study.match_reason = best
        study.match_score = MATCH_SCORES[best[0]]
    else:
        study.match_type, study.match_score, study.match_reason = "none", 0.0, "no exact/synonym trait agreement"
    return study
