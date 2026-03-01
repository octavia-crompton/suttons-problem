#!/usr/bin/env python3
"""
advective_et_litsearch.py

Find studies on advective enhancement of evapotranspiration (ET),
pull abstracts + citation counts, and write:
  - CSV
  - Markdown report

Sources:
  - OpenAlex Works API: cited_by_count; abstract_inverted_index reconstructed
  - Semantic Scholar Graph API: abstract; citationCount; DOI via externalIds.DOI

Notebook-safe:
  - argparse uses parse_known_args() so ipykernel's injected "-f <kernel.json>" won't crash it.

Install:
  pip install requests pandas
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from difflib import SequenceMatcher
from typing import Any, Dict, List, Optional, Tuple

import requests
import pandas as pd


DEFAULT_QUERIES = [
    "advective enhancement evapotranspiration",
    "local advection evapotranspiration",
    "oasis effect evaporation evapotranspiration",
    "fetch effect evaporation",
    "hot dry air advection evaporation",
    "advection irrigation evapotranspiration",
    "internal boundary layer evapotranspiration advection",
    "sensible heat advection latent heat flux",
    "patch scale advection evapotranspiration",
]

KEYWORDS = [
    "advection",
    "advective",
    "local advection",
    "oasis effect",
    "fetch",
    "internal boundary layer",
    "boundary layer",
    "hot dry air",
    "sensible heat",
    "latent heat",
    "heterogeneity",
    "patch",
    "irrigat",
    "evaporative demand",
    "aerodynamic",
    "energy balance",
    "evapotranspiration enhancement",
]


# ---------------------------
# Utilities
# ---------------------------

def _norm_doi(doi: Optional[str]) -> Optional[str]:
    if not doi:
        return None
    doi = doi.strip().lower()
    doi = doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
    doi = doi.replace("doi:", "").strip()
    return doi or None


def _norm_title(t: str) -> str:
    t = t.lower()
    t = re.sub(r"[^a-z0-9\s]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def _sleep_polite(seconds: float):
    if seconds and seconds > 0:
        time.sleep(seconds)


def request_json(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = 30,
    max_retries: int = 5,
    backoff: float = 1.7,
    polite_sleep: float = 0.2,
) -> Dict[str, Any]:
    """
    - Retries on transient failures + 429 rate limiting.
    - DOES NOT retry on permanent 4xx (except 429), and prints the API error body.
    """
    last_err = None
    headers = dict(headers or {})
    headers.setdefault("User-Agent", "advective_et_litsearch/1.1 (requests)")

    for i in range(max_retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)

            if r.status_code == 429:
                wait = (backoff ** i) + 1.0
                _sleep_polite(wait)
                continue

            # Don't retry other 4xx: it's usually a parameter error (like bad fields=)
            if 400 <= r.status_code < 500:
                raise RuntimeError(f"HTTP {r.status_code} for {r.url}\n{r.text[:1500]}")

            r.raise_for_status()
            _sleep_polite(polite_sleep)
            return r.json()

        except Exception as e:
            last_err = e
            wait = (backoff ** i) + 0.5
            _sleep_polite(wait)

    raise RuntimeError(f"Failed request after {max_retries} tries: {url} ({last_err})")


def reconstruct_openalex_abstract(inv_idx: Optional[Dict[str, List[int]]]) -> Optional[str]:
    """
    OpenAlex provides abstract_inverted_index (word -> positions).
    Reconstruct by placing words into their positions and joining with spaces.
    """
    if not inv_idx or not isinstance(inv_idx, dict):
        return None

    max_pos = -1
    for positions in inv_idx.values():
        if positions:
            max_pos = max(max_pos, max(positions))
    if max_pos < 0:
        return None

    tokens = [""] * (max_pos + 1)
    for word, positions in inv_idx.items():
        for p in positions:
            if 0 <= p <= max_pos:
                tokens[p] = word

    text = " ".join(t for t in tokens if t)
    text = re.sub(r"\s+", " ", text).strip()
    return text or None


def relevance_score(title: str, abstract: Optional[str]) -> int:
    txt = ((title or "") + "\n" + (abstract or "")).lower()
    score = 0
    for kw in KEYWORDS:
        if kw in txt:
            score += 2 if " " in kw else 1
    if "advec" in (title or "").lower():
        score += 3
    return score


def summarize_abstract(abstract: Optional[str], max_bullets: int = 3) -> str:
    """
    Lightweight extractive summary: pick high keyword-overlap sentences.
    """
    if not abstract:
        return ""

    sents = re.split(r"(?<=[.!?])\s+", abstract.strip())
    sents = [s.strip() for s in sents if len(s.strip()) >= 30]

    def sent_score(s: str) -> Tuple[int, int]:
        s_l = s.lower()
        k = sum(1 for kw in KEYWORDS if kw in s_l)
        return (k, min(len(s), 240))

    ranked = sorted(sents, key=sent_score, reverse=True)

    chosen: List[str] = []
    for s in ranked:
        if len(chosen) >= max_bullets:
            break
        if any(SequenceMatcher(None, s, c).ratio() > 0.85 for c in chosen):
            continue
        chosen.append(s)

    if not chosen and sents:
        chosen = sents[:max_bullets]

    return "\n".join(f"- {c}" for c in chosen)


# ---------------------------
# OpenAlex
# ---------------------------

def openalex_search(
    query: str,
    max_results: int,
    email: Optional[str],
    polite_sleep: float,
) -> List[Dict[str, Any]]:
    base = "https://api.openalex.org/works"
    params = {
        "search": query,
        "per-page": 200,
        "sort": "cited_by_count:desc",
    }
    if email:
        params["mailto"] = email

    out: List[Dict[str, Any]] = []
    cursor = "*"

    while len(out) < max_results:
        params2 = dict(params)
        params2["cursor"] = cursor
        data = request_json(base, params=params2, polite_sleep=polite_sleep)
        results = data.get("results") or []
        out.extend(results)

        cursor = data.get("meta", {}).get("next_cursor")
        if not cursor or not results:
            break

    return out[:max_results]


def openalex_to_record(w: Dict[str, Any], query_tag: str) -> Dict[str, Any]:
    doi = _norm_doi(w.get("doi"))
    title = (w.get("title") or "").strip()
    year = w.get("publication_year")
    cited_by = w.get("cited_by_count")

    venue = None
    try:
        venue = w.get("primary_location", {}).get("source", {}).get("display_name")
    except Exception:
        venue = None

    authors = []
    for a in (w.get("authorships") or [])[:12]:
        nm = a.get("author", {}).get("display_name")
        if nm:
            authors.append(nm)
    authors_str = ", ".join(authors) if authors else None

    abstract = reconstruct_openalex_abstract(w.get("abstract_inverted_index"))
    url = w.get("id")  # OpenAlex work URL

    return {
        "title": title or None,
        "year": year,
        "doi": doi,
        "venue": venue,
        "authors": authors_str,
        "abstract": abstract,
        "abstract_source": "openalex" if abstract else None,
        "openalex_cited_by_count": cited_by,
        "s2_citation_count": None,
        "url": url,
        "query_hits": query_tag,
    }


# ---------------------------
# Semantic Scholar
# ---------------------------

def s2_search(
    query: str,
    max_results: int,
    api_key: Optional[str],
    polite_sleep: float,
) -> List[Dict[str, Any]]:
    """
    Semantic Scholar search endpoint.
    Paginates using offset/limit if max_results > 100.
    """
    base = "https://api.semanticscholar.org/graph/v1/paper/search"
    headers: Dict[str, str] = {}
    if api_key:
        headers["x-api-key"] = api_key

    # IMPORTANT: DOI is in externalIds.DOI (not a top-level "doi" field for this endpoint)
    fields = "title,year,abstract,citationCount,externalIds,url,venue,authors"

    out: List[Dict[str, Any]] = []
    offset = 0

    while len(out) < max_results:
        limit = min(100, max_results - len(out))
        params = {
            "query": query,
            "limit": limit,
            "offset": offset,
            "fields": fields,
        }
        data = request_json(base, params=params, headers=headers, polite_sleep=polite_sleep)
        batch = data.get("data") or []
        if not batch:
            break
        out.extend(batch)
        offset += len(batch)

        # If API returns fewer than requested, we likely hit the end
        if len(batch) < limit:
            break

    return out[:max_results]


def s2_to_record(p: Dict[str, Any], query_tag: str) -> Dict[str, Any]:
    ext = p.get("externalIds") or {}
    doi = _norm_doi(ext.get("DOI") or ext.get("doi"))

    title = (p.get("title") or "").strip()
    year = p.get("year")
    cited = p.get("citationCount")
    abstract = (p.get("abstract") or "").strip() or None
    venue = p.get("venue")
    url = p.get("url")

    authors = []
    for a in (p.get("authors") or [])[:12]:
        nm = a.get("name")
        if nm:
            authors.append(nm)
    authors_str = ", ".join(authors) if authors else None

    return {
        "title": title or None,
        "year": year,
        "doi": doi,
        "venue": venue,
        "authors": authors_str,
        "abstract": abstract,
        "abstract_source": "semanticscholar" if abstract else None,
        "openalex_cited_by_count": None,
        "s2_citation_count": cited,
        "url": url,
        "query_hits": query_tag,
    }


# ---------------------------
# Merge logic
# ---------------------------

def merge_records(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge two records:
      - Prefer Semantic Scholar for abstract and URL if present
      - Keep both citation counts
      - Union query_hits
    """
    out = dict(a)

    # Prefer b for empty fields
    for k, v in b.items():
        if k in ("openalex_cited_by_count", "s2_citation_count"):
            continue
        if out.get(k) in (None, "", []) and v not in (None, "", []):
            out[k] = v

    # citation counts
    if out.get("openalex_cited_by_count") is None and b.get("openalex_cited_by_count") is not None:
        out["openalex_cited_by_count"] = b["openalex_cited_by_count"]
    if out.get("s2_citation_count") is None and b.get("s2_citation_count") is not None:
        out["s2_citation_count"] = b["s2_citation_count"]

    # abstract preference: Semantic Scholar > OpenAlex
    if b.get("abstract_source") == "semanticscholar" and b.get("abstract"):
        out["abstract"] = b["abstract"]
        out["abstract_source"] = "semanticscholar"

    # URL preference
    if b.get("url"):
        out["url"] = b["url"]

    # query hits: concatenate unique
    ah = set((out.get("query_hits") or "").split(" | ")) if out.get("query_hits") else set()
    bh = set((b.get("query_hits") or "").split(" | ")) if b.get("query_hits") else set()
    hits = sorted([h for h in (ah | bh) if h])
    out["query_hits"] = " | ".join(hits) if hits else None

    return out


def fuzzy_join_no_doi(
    openalex: List[Dict[str, Any]],
    s2: List[Dict[str, Any]],
    title_threshold: float = 0.92,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    For records without DOI, attempt title+year fuzzy match.
    Returns merged list + leftovers list.
    """
    s2_pool = [r for r in s2 if not r.get("doi") and r.get("title")]
    used = set()

    merged: List[Dict[str, Any]] = []
    leftovers: List[Dict[str, Any]] = []

    for oa in openalex:
        if oa.get("doi") or not oa.get("title"):
            leftovers.append(oa)
            continue

        best_j = None
        best_score = 0.0
        oa_t = _norm_title(oa["title"])
        oa_y = oa.get("year")

        for j, srec in enumerate(s2_pool):
            if j in used:
                continue
            if oa_y and srec.get("year") and abs(int(oa_y) - int(srec["year"])) > 1:
                continue
            sc = SequenceMatcher(None, oa_t, _norm_title(srec["title"])).ratio()
            if sc > best_score:
                best_score = sc
                best_j = j

        if best_j is not None and best_score >= title_threshold:
            used.add(best_j)
            merged.append(merge_records(oa, s2_pool[best_j]))
        else:
            leftovers.append(oa)

    for j, srec in enumerate(s2_pool):
        if j not in used:
            leftovers.append(srec)

    return merged, leftovers


# ---------------------------
# Output
# ---------------------------

def write_markdown(df: pd.DataFrame, out_md: str, top_n: int = 50):
    lines: List[str] = []
    lines.append("# Advective enhancement of ET — literature scan")
    lines.append("")
    lines.append(f"Records: {len(df)}")
    lines.append("")

    show = df.head(top_n).copy()

    for _, row in show.iterrows():
        title = row.get("title") or "(no title)"
        year = row.get("year") or ""
        doi = row.get("doi") or ""
        venue = row.get("venue") or ""
        authors = row.get("authors") or ""
        url = row.get("url") or ""
        c_oa = row.get("openalex_cited_by_count")
        c_s2 = row.get("s2_citation_count")

        lines.append(f"## {title}")
        meta = []
        if year:
            meta.append(str(year))
        if venue:
            meta.append(str(venue))
        if authors:
            meta.append(str(authors))
        if meta:
            lines.append("*" + " — ".join(meta) + "*")
        if doi:
            lines.append(f"- DOI: {doi}")
        if url:
            lines.append(f"- Link: {url}")
        lines.append(f"- Citations (OpenAlex): {c_oa if pd.notna(c_oa) else 'NA'}")
        lines.append(f"- Citations (Semantic Scholar): {c_s2 if pd.notna(c_s2) else 'NA'}")

        summ = row.get("abstract_summary") or ""
        if summ:
            lines.append("")
            lines.append("**Abstract summary (extractive):**")
            lines.append(summ)
        lines.append("")

    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# ---------------------------
# Core pipeline (callable from notebooks)
# ---------------------------

def run(
    email: Optional[str] = None,
    s2_api_key: Optional[str] = None,
    max_per_query: int = 80,
    min_year: Optional[int] = None,
    max_year: Optional[int] = None,
    min_relevance: int = 2,
    polite_sleep: float = 0.2,
    out_csv: str = "advective_et_literature.csv",
    out_md: str = "advective_et_literature.md",
    queries: Optional[List[str]] = None,
    queries_json: Optional[str] = None,
    save: bool = True,
    md_top_n: int = 50,
    use_semantic_scholar: bool = True,
) -> pd.DataFrame:
    """
    Notebook-friendly entry point. Returns a DataFrame.
    If save=True, also writes CSV + MD to disk.
    """
    qlist = list(queries) if queries else list(DEFAULT_QUERIES)
    if queries_json:
        with open(queries_json, "r", encoding="utf-8") as f:
            obj = json.load(f)
        qlist = obj.get("queries") or qlist

    s2_key = s2_api_key or os.environ.get("S2_API_KEY")

    openalex_records: List[Dict[str, Any]] = []
    s2_records: List[Dict[str, Any]] = []

    for q in qlist:
        oa = openalex_search(q, max_results=max_per_query, email=email, polite_sleep=polite_sleep)
        openalex_records.extend(openalex_to_record(w, q) for w in oa)

        if use_semantic_scholar:
            try:
                s2 = s2_search(q, max_results=max_per_query, api_key=s2_key, polite_sleep=polite_sleep)
                s2_records.extend(s2_to_record(p, q) for p in s2)
            except Exception as e:
                print(f"[warn] Semantic Scholar failed for query={q!r}: {e}")

    # Deduplicate within each source
    def dedupe(recs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        by_doi: Dict[str, Dict[str, Any]] = {}
        by_key: Dict[Tuple[str, Any], Dict[str, Any]] = {}
        for r in recs:
            doi = r.get("doi")
            if doi:
                by_doi[doi] = merge_records(by_doi[doi], r) if doi in by_doi else r
            else:
                key = (_norm_title(r.get("title") or ""), r.get("year"))
                if key[0]:
                    by_key[key] = merge_records(by_key[key], r) if key in by_key else r
        return list(by_doi.values()) + list(by_key.values())

    openalex_u = dedupe(openalex_records)
    s2_u = dedupe(s2_records) if use_semantic_scholar else []

    # Merge across sources by DOI
    merged_by_doi: Dict[str, Dict[str, Any]] = {}
    no_doi_oa: List[Dict[str, Any]] = []
    no_doi_s2: List[Dict[str, Any]] = []

    for r in openalex_u:
        if r.get("doi"):
            merged_by_doi[r["doi"]] = merge_records(merged_by_doi.get(r["doi"], r), r)
        else:
            no_doi_oa.append(r)

    for r in s2_u:
        if r.get("doi"):
            merged_by_doi[r["doi"]] = merge_records(merged_by_doi.get(r["doi"], r), r)
        else:
            no_doi_s2.append(r)

    merged: List[Dict[str, Any]] = list(merged_by_doi.values())

    # Fuzzy merge remaining no-DOI
    if use_semantic_scholar:
        fuzzy_merged, leftovers = fuzzy_join_no_doi(no_doi_oa, no_doi_s2)
        merged.extend(fuzzy_merged)
        merged.extend(leftovers)
    else:
        merged.extend(no_doi_oa)

    # Year filters
    if min_year is not None:
        merged = [r for r in merged if (r.get("year") is None or int(r["year"]) >= min_year)]
    if max_year is not None:
        merged = [r for r in merged if (r.get("year") is None or int(r["year"]) <= max_year)]

    # Relevance + summaries
    for r in merged:
        r["relevance"] = relevance_score(r.get("title") or "", r.get("abstract"))
        r["abstract_summary"] = summarize_abstract(r.get("abstract"))

    merged = [r for r in merged if r.get("relevance", 0) >= min_relevance]

    # Sort: relevance desc, then max citations desc
    def max_cites(r: Dict[str, Any]) -> int:
        c1 = r.get("openalex_cited_by_count") or 0
        c2 = r.get("s2_citation_count") or 0
        try:
            return int(max(c1, c2))
        except Exception:
            return 0

    merged.sort(key=lambda r: (r.get("relevance") or 0, max_cites(r)), reverse=True)

    df = pd.DataFrame(merged)

    cols = [
        "title", "year", "authors", "venue", "doi", "url",
        "openalex_cited_by_count", "s2_citation_count",
        "relevance", "abstract_source", "query_hits",
        "abstract_summary", "abstract",
    ]
    for c in cols:
        if c not in df.columns:
            df[c] = None
    df = df[cols]

    if save:
        df.to_csv(out_csv, index=False)
        write_markdown(df, out_md, top_n=md_top_n)

    return df


# ---------------------------
# CLI (still works; notebook-safe parse)
# ---------------------------

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", default=None, help="Email for OpenAlex 'mailto=' parameter")
    ap.add_argument("--s2-api-key", default=None, help="Semantic Scholar API key (or set S2_API_KEY env var)")
    ap.add_argument("--max-per-query", type=int, default=80, help="Max results per query per source")
    ap.add_argument("--min-year", type=int, default=None)
    ap.add_argument("--max-year", type=int, default=None)
    ap.add_argument("--min-relevance", type=int, default=2)
    ap.add_argument("--polite-sleep", type=float, default=0.2)
    ap.add_argument("--out-csv", default="advective_et_literature.csv")
    ap.add_argument("--out-md", default="advective_et_literature.md")
    ap.add_argument("--queries-json", default=None, help="JSON file: {\"queries\": [\"...\"]}")
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--no-s2", action="store_true", help="Disable Semantic Scholar (OpenAlex only)")

    # Notebook fix: ignores ipykernel's injected args like "-f <kernel.json>"
    args, _unknown = ap.parse_known_args(argv)

    df = run(
        email=args.email,
        s2_api_key=args.s2_api_key,
        max_per_query=args.max_per_query,
        min_year=args.min_year,
        max_year=args.max_year,
        min_relevance=args.min_relevance,
        polite_sleep=args.polite_sleep,
        out_csv=args.out_csv,
        out_md=args.out_md,
        queries_json=args.queries_json,
        save=(not args.no_save),
        use_semantic_scholar=(not args.no_s2),
    )

    print(f"Rows: {len(df)}")
    if not args.no_save:
        print(f"Wrote: {args.out_csv}")
        print(f"Wrote: {args.out_md}")


if __name__ == "__main__":
    main()
