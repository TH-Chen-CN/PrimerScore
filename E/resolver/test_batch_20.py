"""
Batch smoke test for resolver using D module output.

Usage:
  python -m E.resolver.test_batch_20 --d_output D_output.json
Optional:
  --n 20
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, List

from E.resolver.resolver import TranscriptResolver, strip_version


def collect_subject_ids(d_results: List[Dict[str, Any]]) -> List[str]:
    ids: List[str] = []
    for pair in d_results:
        tid = pair.get("target_subject_id")
        if tid:
            ids.append(str(tid))

        off = (pair.get("off_targets") or {}).get("near") or []
        for amp in off:
            sid = amp.get("subject_id")
            if sid:
                ids.append(str(sid))

    # dedupe keep order + strip version
    seen = set()
    uniq: List[str] = []
    for x in ids:
        x0 = strip_version(x)
        if x0 and x0 not in seen:
            seen.add(x0)
            uniq.append(x0)
    return uniq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d_output", required=True, help="Path to D output JSON file (list[dict])")
    ap.add_argument("--n", type=int, default=20, help="How many unique subject_ids to test")
    ap.add_argument("--cache", default="resolver_cache.sqlite", help="Cache sqlite path")
    ap.add_argument("--no_ensembl", action="store_true", help="Disable Ensembl fallback")
    args = ap.parse_args()

    with open(args.d_output, "r", encoding="utf-8") as f:
        d_results = json.load(f)

    if not isinstance(d_results, list):
        raise ValueError("D output must be a JSON list[dict].")

    ids = collect_subject_ids(d_results)
    test_ids = ids[: args.n]

    print("Testing IDs:", test_ids)

    resolver = TranscriptResolver(cache_path=args.cache, enable_ensembl_fallback=(not args.no_ensembl))
    res = resolver.resolve_many(test_ids)

    # Stats
    total = len(res)
    hit = sum(1 for r in res if r.get("status") == "hit" and r.get("glyma_id"))
    amb = sum(1 for r in res if r.get("status") == "ambiguous" and r.get("glyma_id"))
    miss = sum(1 for r in res if r.get("status") == "miss")
    err = sum(1 for r in res if r.get("status") == "error")

    print(f"\nSummary: total={total} hit={hit} ambiguous={amb} miss={miss} error={err}\n")

    # Print results (compact)
    for r in res:
        tid = r.get("transcript_id")
        gid = r.get("glyma_id")
        st = r.get("status")
        src = r.get("source")
        gene_id = r.get("gene_id")
        locus = r.get("raw_locus_tag")
        notes = r.get("notes")
        print(f"- {tid:>15}  status={st:<10}  glyma={gid}  source={src}  gene_id={gene_id}  locus={locus}  notes={notes}")

    # Miss list for debugging
    print("\n[Miss/Error items]")
    for r in res:
        if not r.get("glyma_id") or r.get("status") in ("miss", "error"):
            print(r)


if __name__ == "__main__":
    main()
