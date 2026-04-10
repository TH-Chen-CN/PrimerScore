"""
Resolver: XM_/NM_ (RefSeq transcript accession) -> Glyma gene model ID

Primary strategy (A):
- NCBI E-utilities:
  XM/NM accession -> esearch nuccore UID -> elink gene -> esummary gene
  Extract locus tag like "GLYMA_09G234900v4" from gene record
  Normalize to "Glyma.09G234900"

Fallback (optional, B):
- Ensembl Plants REST xrefs (best-effort)

Requirements:
- Batch-friendly
- SQLite cache
- Explainable output fields (gene_id, raw_locus_tag, etc.)

Run tests:
  python -m E.resolver.test_batch_20 --d_output D_output.json
"""

from __future__ import annotations

import json
import re
import sqlite3
import time
from dataclasses import dataclass
from typing import Dict, List, Optional

import requests


# =========================
# Config
# =========================

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_BASE = "https://rest.ensembl.org"

DEFAULT_DB_PATH = "resolver_cache.sqlite"
DEFAULT_TIMEOUT = 15
USER_AGENT = "primer-e-resolver/0.2"

# Minimal network politeness (NCBI can throttle)
DEFAULT_SLEEP_BETWEEN_CALLS = 0.15


# =========================
# Utilities
# =========================

def strip_version(accession: str) -> str:
    """NM_001234.2 -> NM_001234"""
    return accession.split(".")[0].strip()


# Keep this (sometimes GenBank text contains Glyma.XXGXXXXXX directly)
def extract_glyma(text: str) -> List[str]:
    """Extract Glyma IDs from free text. Matches: Glyma.12G040000"""
    if not text:
        return []
    return list(set(re.findall(r"Glyma\.\d{2}G\d{6}", text)))


# NEW: locus tag extraction from NCBI Gene record (your screenshot shows this!)
GLYMA_LOCUS_RE = re.compile(r"\bGLYMA_\d{2}G\d{6}(?:v\d+)?\b")


def extract_glyma_locus(text: str) -> List[str]:
    """Extract NCBI locus tag like GLYMA_09G234900v4 or GLYMA_09G234900."""
    if not text:
        return []
    return list(set(GLYMA_LOCUS_RE.findall(text)))


def normalize_ncbi_locus_to_glyma(locus: str) -> Optional[str]:
    """
    GLYMA_09G234900v4 -> Glyma.09G234900
    GLYMA_09G234900   -> Glyma.09G234900
    """
    if not locus:
        return None
    x = locus.strip()
    x = re.sub(r"v\d+$", "", x)              # drop version suffix
    x = x.replace("GLYMA_", "Glyma_", 1)     # prefix case
    x = x.replace("Glyma_", "Glyma.", 1)     # delimiter
    if re.fullmatch(r"Glyma\.\d{2}G\d{6}", x):
        return x
    return None


def _http_get(url: str, timeout: int = DEFAULT_TIMEOUT) -> requests.Response:
    """Simple GET with UA + basic throttle."""
    time.sleep(DEFAULT_SLEEP_BETWEEN_CALLS)
    return requests.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})


# =========================
# Cache
# =========================

class ResolverCache:
    """
    SQLite cache for resolver results.

    Schema is auto-migrated if you already have an old resolver_cache.sqlite.
    """
    def __init__(self, db_path: str = DEFAULT_DB_PATH):
        self.db_path = db_path
        self._init_db_and_migrate()

    def _init_db_and_migrate(self):
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()

        # Create table if not exists (new schema)
        cur.execute("""
        CREATE TABLE IF NOT EXISTS resolver_cache (
            transcript_id TEXT PRIMARY KEY,
            glyma_id TEXT,
            status TEXT,
            source TEXT,
            gene_id TEXT,
            raw_locus_tag TEXT,
            timestamp REAL
        )
        """)
        conn.commit()

        # Migration: add missing columns if older schema exists
        cur.execute("PRAGMA table_info(resolver_cache)")
        cols = {row[1] for row in cur.fetchall()}

        if "gene_id" not in cols:
            cur.execute("ALTER TABLE resolver_cache ADD COLUMN gene_id TEXT")
        if "raw_locus_tag" not in cols:
            cur.execute("ALTER TABLE resolver_cache ADD COLUMN raw_locus_tag TEXT")
        if "timestamp" not in cols:
            cur.execute("ALTER TABLE resolver_cache ADD COLUMN timestamp REAL")
        conn.commit()
        conn.close()

    def get(self, transcript_id: str) -> Optional[Dict]:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute(
            "SELECT glyma_id, status, source, gene_id, raw_locus_tag FROM resolver_cache WHERE transcript_id=?",
            (transcript_id,)
        )
        row = cur.fetchone()
        conn.close()
        if not row:
            return None
        return {
            "glyma_id": row[0],
            "status": row[1],
            "source": row[2],
            "gene_id": row[3],
            "raw_locus_tag": row[4],
        }

    def set(
        self,
        transcript_id: str,
        glyma_id: Optional[str],
        status: str,
        source: str,
        gene_id: Optional[str] = None,
        raw_locus_tag: Optional[str] = None,
    ) -> None:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute("""
        INSERT OR REPLACE INTO resolver_cache
        (transcript_id, glyma_id, status, source, gene_id, raw_locus_tag, timestamp)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (transcript_id, glyma_id, status, source, gene_id, raw_locus_tag, time.time()))
        conn.commit()
        conn.close()


# =========================
# Resolver result structure
# =========================

@dataclass
class ResolveResult:
    transcript_id: str
    glyma_id: Optional[str]
    status: str                 # hit | ambiguous | miss | error
    source: str                 # NCBI_GENE | NCBI_GB | ENSEMBL
    gene_id: Optional[str] = None
    raw_locus_tag: Optional[str] = None
    notes: Optional[str] = None

    def to_dict(self) -> Dict:
        return {
            "transcript_id": self.transcript_id,
            "glyma_id": self.glyma_id,
            "status": self.status,
            "source": self.source,
            "gene_id": self.gene_id,
            "raw_locus_tag": self.raw_locus_tag,
            "notes": self.notes,
        }


# =========================
# NCBI Resolver (A)
# =========================

def ncbi_resolve(accession: str) -> ResolveResult:
    """
    Primary: NCBI Gene route.
    accession -> esearch nuccore UID -> elink gene -> esummary gene -> extract GLYMA locus tag -> normalize

    (Optional) We keep a "GB text scan" pre-step as a quick win:
    accession -> efetch gb -> scan Glyma.XXGXXXXXX in features
    """
    acc = strip_version(accession)

    try:
        # 0) esearch: accession -> nuccore UID
        search_url = f"{NCBI_BASE}/esearch.fcgi?db=nuccore&term={acc}[Accession]&retmode=json"
        r = _http_get(search_url)
        r.raise_for_status()
        j = r.json()
        idlist = j.get("esearchresult", {}).get("idlist", [])

        if not idlist:
            # fallback: broader search
            search_url2 = f"{NCBI_BASE}/esearch.fcgi?db=nuccore&term={acc}&retmode=json"
            r = _http_get(search_url2)
            r.raise_for_status()
            j = r.json()
            idlist = j.get("esearchresult", {}).get("idlist", [])

        if not idlist:
            return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="NCBI_GENE",
                                 notes="esearch nuccore idlist empty")

        uid = idlist[0]

        # 1) quick win: efetch GenBank text and scan Glyma.XXGXXXXXX directly (not always present)
        fetch_url = f"{NCBI_BASE}/efetch.fcgi?db=nuccore&id={uid}&rettype=gb&retmode=text"
        r2 = _http_get(fetch_url)
        r2.raise_for_status()
        gb_text = r2.text

        glymas = extract_glyma(gb_text)
        if len(glymas) == 1:
            return ResolveResult(transcript_id=accession, glyma_id=glymas[0], status="hit", source="NCBI_GB",
                                 notes="found Glyma.* in GenBank text")
        elif len(glymas) > 1:
            return ResolveResult(transcript_id=accession, glyma_id=sorted(glymas)[0], status="ambiguous", source="NCBI_GB",
                                 notes=f"multiple Glyma.* in GenBank text: {sorted(glymas)[:5]}")

        # 2) main route: nuccore UID -> gene link -> gene summary -> parse GLYMA locus tag
        link_url = f"{NCBI_BASE}/elink.fcgi?dbfrom=nuccore&db=gene&id={uid}&retmode=json"
        r3 = _http_get(link_url)
        r3.raise_for_status()
        data = r3.json()

        gene_ids: List[str] = []
        for ls in data.get("linksets", []):
            for ldb in ls.get("linksetdbs", []):
                if ldb.get("dbto") == "gene":
                    gene_ids.extend([str(x) for x in ldb.get("links", [])])

        if not gene_ids:
            return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="NCBI_GENE",
                                 notes="elink nuccore->gene returned empty")

        gid = gene_ids[0]

        sum_url = f"{NCBI_BASE}/esummary.fcgi?db=gene&id={gid}&retmode=json"
        r4 = _http_get(sum_url)
        r4.raise_for_status()
        sdata = r4.json()

        rec = sdata.get("result", {}).get(str(gid), {})
        text_blob = json.dumps(rec)

        loci = extract_glyma_locus(text_blob)
        if not loci:
            return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="NCBI_GENE",
                                 gene_id=str(gid), notes="no GLYMA locus tag found in gene record")

        loci_sorted = sorted(loci)
        raw_locus = loci_sorted[0]
        glyma_norm = normalize_ncbi_locus_to_glyma(raw_locus)

        if not glyma_norm:
            return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="NCBI_GENE",
                                 gene_id=str(gid), raw_locus_tag=raw_locus,
                                 notes="found GLYMA locus tag but failed to normalize")

        if len(loci_sorted) == 1:
            return ResolveResult(transcript_id=accession, glyma_id=glyma_norm, status="hit", source="NCBI_GENE",
                                 gene_id=str(gid), raw_locus_tag=raw_locus,
                                 notes="normalized from NCBI Gene locus tag")
        else:
            return ResolveResult(transcript_id=accession, glyma_id=glyma_norm, status="ambiguous", source="NCBI_GENE",
                                 gene_id=str(gid), raw_locus_tag=raw_locus,
                                 notes=f"multiple locus tags found: {loci_sorted[:5]}")

    except Exception as e:
        return ResolveResult(transcript_id=accession, glyma_id=None, status="error", source="NCBI_GENE",
                             notes=f"exception: {type(e).__name__}: {e}")


# =========================
# Ensembl Resolver (B) - best effort fallback
# =========================

def ensembl_resolve(accession: str) -> ResolveResult:
    """
    Best-effort fallback: Ensembl REST xrefs.
    Might miss depending on Ensembl data coverage.
    """
    acc = strip_version(accession)
    try:
        headers = {
            "Content-Type": "application/json",
            "User-Agent": USER_AGENT,
        }
        url = f"{ENSEMBL_BASE}/xrefs/name/glycine_max/{acc}"
        r = requests.get(url, headers=headers, timeout=DEFAULT_TIMEOUT)
        if r.status_code != 200:
            return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="ENSEMBL",
                                 notes=f"http {r.status_code}")

        hits = r.json()
        glymas: List[str] = []
        for h in hits:
            for v in h.values():
                if isinstance(v, str):
                    glymas.extend(extract_glyma(v))
        glymas = sorted(set(glymas))

        if len(glymas) == 1:
            return ResolveResult(transcript_id=accession, glyma_id=glymas[0], status="hit", source="ENSEMBL")
        if len(glymas) > 1:
            return ResolveResult(transcript_id=accession, glyma_id=glymas[0], status="ambiguous", source="ENSEMBL",
                                 notes=f"multiple glyma hits: {glymas[:5]}")
        return ResolveResult(transcript_id=accession, glyma_id=None, status="miss", source="ENSEMBL")

    except Exception as e:
        return ResolveResult(transcript_id=accession, glyma_id=None, status="error", source="ENSEMBL",
                             notes=f"exception: {type(e).__name__}: {e}")


# =========================
# Public API
# =========================

class TranscriptResolver:
    """
    Public resolver object:
      - resolve_one()
      - resolve_many()
    """

    def __init__(self, cache_path: str = DEFAULT_DB_PATH, enable_ensembl_fallback: bool = True):
        self.cache = ResolverCache(cache_path)
        self.enable_ensembl_fallback = enable_ensembl_fallback

    def resolve_one(self, transcript_id: str) -> Dict:
        """
        Resolve a single XM_/NM_ accession.
        Returns a dict with explainable fields.
        """
        tid = strip_version(transcript_id)

        # Cache lookup (cache by stripped version)
        cached = self.cache.get(tid)
        if cached:
            return ResolveResult(
                transcript_id=tid,
                glyma_id=cached.get("glyma_id"),
                status=cached.get("status") or "miss",
                source=cached.get("source") or "CACHE",
                gene_id=cached.get("gene_id"),
                raw_locus_tag=cached.get("raw_locus_tag"),
                notes="cache_hit",
            ).to_dict()

        # A: NCBI
        res = ncbi_resolve(tid)

        # B: Ensembl fallback (only if miss/error OR no glyma_id)
        if self.enable_ensembl_fallback and (not res.glyma_id or res.status in ("miss", "error")):
            res2 = ensembl_resolve(tid)
            # Use fallback only if it actually finds a glyma_id
            if res2.glyma_id:
                res = res2

        # Write cache (cache final result)
        self.cache.set(
            transcript_id=tid,
            glyma_id=res.glyma_id,
            status=res.status,
            source=res.source,
            gene_id=res.gene_id,
            raw_locus_tag=res.raw_locus_tag,
        )
        return res.to_dict()

    def resolve_many(self, transcript_ids: List[str]) -> List[Dict]:
        """
        Batch resolve with de-dup (keep order).
        """
        seen = set()
        uniq: List[str] = []
        for x in transcript_ids:
            x0 = strip_version(x)
            if x0 and x0 not in seen:
                seen.add(x0)
                uniq.append(x0)

        out: List[Dict] = []
        for tid in uniq:
            out.append(self.resolve_one(tid))
        return out
