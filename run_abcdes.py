from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from pipeline_ABCDES import pipeline_ABCDES


# 你自己的目标序列放这里
SEQUENCE = """
ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTGTTACACAAGGTGAATCTTGGAGTGATTTTCTTCTTGACCTACAAGGCTAG
"""

OUTPUT_JSON = "abcdes_result.json"


def save_json(path: str | Path, obj: Any) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def print_meta_summary(result: Dict[str, Any]) -> None:
    meta = result.get("meta", {}) or {}
    counts = meta.get("counts", {}) or {}
    e_status = meta.get("e_status", {}) or {}
    scoring_status = meta.get("scoring_status", {}) or {}

    print("\n========== PIPELINE SUMMARY ==========")
    print(f"Pipeline: {meta.get('pipeline')}")
    print(f"Target length: {meta.get('target_length')}")
    print(f"Organism: {meta.get('organism')}")
    print(f"DB: {meta.get('db')}")

    print("\n[Counts]")
    for k, v in counts.items():
        print(f"  {k}: {v}")

    print("\n[E status]")
    print(f"  enabled: {e_status.get('enabled')}")
    print(f"  ok: {e_status.get('ok')}")
    print(f"  error: {e_status.get('error')}")
    print(f"  resolver_enabled: {e_status.get('resolver_enabled')}")
    print(f"  resolved_target_pairs: {e_status.get('resolved_target_pairs')}")
    print(f"  pairs_with_nonempty_target_expr: {e_status.get('pairs_with_nonempty_target_expr')}")
    print(f"  expression_sources_used: {e_status.get('expression_sources_used')}")

    print("\n[Scoring status]")
    print(f"  enabled: {scoring_status.get('enabled')}")
    print(f"  ok: {scoring_status.get('ok')}")
    print(f"  error: {scoring_status.get('error')}")
    print(f"  n_pairs_in_d: {scoring_status.get('n_pairs_in_d')}")
    print(f"  n_pairs_in_e: {scoring_status.get('n_pairs_in_e')}")


def print_top_thermo(rankings: Dict[str, Any], top_n: int = 5) -> None:
    rows: List[Dict[str, Any]] = rankings.get("thermo_ranking", []) or []

    print("\n========== TOP THERMO ==========")
    if not rows:
        print("No thermo ranking results.")
        return

    for i, row in enumerate(rows[:top_n], 1):
        print(
            f"{i}. "
            f"pair={row.get('pair_index')} | "
            f"thermo_score={row.get('thermo_score')} | "
            f"thermo_label={row.get('thermo_label')} | "
            f"specificity_label={row.get('specificity_label')}"
        )
        warnings = row.get("warnings", []) or []
        if warnings:
            for w in warnings:
                print(f"    - {w}")


def print_top_expression(rankings: Dict[str, Any], top_n: int = 10) -> None:
    rows: List[Dict[str, Any]] = rankings.get("expression_ranking", []) or []

    print("\n========== TOP EXPRESSION ==========")
    if not rows:
        print("No expression ranking results.")
        return

    for i, row in enumerate(rows[:top_n], 1):
        print(
            f"{i}. "
            f"pair={row.get('pair_index')} | "
            f"source={row.get('source')} | "
            f"expression_score={row.get('expression_score')} | "
            f"label={row.get('recommendation_label')} | "
            f"target_expr={row.get('target_expr')} | "
            f"near_expr={row.get('near_expr')} | "
            f"far_expr={row.get('far_expr')} | "
            f"target/near={row.get('target_vs_near_ratio')} | "
            f"target/far={row.get('target_vs_far_ratio')}"
        )

        tags = row.get("source_tags", []) or []
        warnings = row.get("source_warnings", []) or []

        if tags:
            print(f"    tags: {', '.join(tags)}")
        if warnings:
            for w in warnings:
                print(f"    - {w}")


def print_top_combined(rankings: Dict[str, Any], top_n: int = 10) -> None:
    rows: List[Dict[str, Any]] = rankings.get("combined_ranking", []) or []

    print("\n========== TOP COMBINED ==========")
    if not rows:
        print("No combined ranking results.")
        return

    for i, row in enumerate(rows[:top_n], 1):
        print(
            f"{i}. "
            f"pair={row.get('pair_index')} | "
            f"source={row.get('source')} | "
            f"combined_score={row.get('combined_score')} | "
            f"expression_score={row.get('expression_score')} | "
            f"thermo_score={row.get('thermo_score')} | "
            f"final={row.get('final_recommendation')} | "
            f"thermo_label={row.get('thermo_label')} | "
            f"specificity_label={row.get('specificity_label')}"
        )
        print(
            f"    target_expr={row.get('target_expr')}, "
            f"near_expr={row.get('near_expr')}, "
            f"far_expr={row.get('far_expr')}, "
            f"target/near={row.get('target_vs_near_ratio')}, "
            f"target/far={row.get('target_vs_far_ratio')}"
        )

        tags = row.get("all_tags", []) or []
        warnings = row.get("warnings", []) or []

        if tags:
            print(f"    tags: {', '.join(tags)}")
        if warnings:
            for w in warnings:
                print(f"    - {w}")


def print_best_recommendations(rankings: Dict[str, Any], top_n: int = 5) -> None:
    rows: List[Dict[str, Any]] = rankings.get("combined_ranking", []) or []

    print("\n========== FINAL RECOMMENDATIONS ==========")
    if not rows:
        print("No final recommendations available.")
        return

    for i, row in enumerate(rows[:top_n], 1):
        print(
            f"{i}. Recommend Pair #{row.get('pair_index')} "
            f"for source/tissue [{row.get('source')}] "
            f"with combined_score={row.get('combined_score')} "
            f"({row.get('final_recommendation')})"
        )


def main() -> None:
    sequence = SEQUENCE.strip()

    result = pipeline_ABCDES(
        sequence=sequence,

        # A
        size_range=(18, 25),

        # B
        min_tm=50.0,
        max_tm=65.0,
        max_hairpin_dG=-3.0,
        max_dimer_dG=-9.0,
        min_gc=0.0,
        max_gc=1.0,
        min_len=18,
        max_len=25,

        # C
        target_length=None,
        max_tm_diff=2.0,
        max_heterodimer_dG=0.0,

        # D (联网)
        organism="Glycine max",
        db="refseq_rna",
        max_workers=3,
        sleep_between_requests=0.5,
        require_strict_orientation=True,
        min_amplicon_len=50,
        max_amplicon_len=3000,
        near_bp=30,
        target_subject_id=None,
        dummy_blast=False,
        min_identity_min=None,

        # E
        enable_e=True,
        expression_cache_db_path="expression_cache.sqlite",
        expression_ttl_seconds=7 * 24 * 3600,
        allowed_canon_conditions=None,
        conekt_prefix="glyma.Wm82.gnm2.ann1.",
        conekt_sleep_s=0.2,
        conekt_timeout=25,
        e_use_cache=True,
        preferred_expression_source="integrated",
        enable_soyomics_fallback=True,
        soyomics_timeout=25,
        source_weights=None,

        # resolver
        resolver_cache_db_path="resolver_cache.sqlite",
        enable_ensembl_fallback=True,

        # S
        enable_scoring=True,
    )

    save_json(OUTPUT_JSON, result)

    print_meta_summary(result)

    rankings = result.get("rankings", {}) or {}
    print_top_thermo(rankings, top_n=5)
    print_top_expression(rankings, top_n=10)
    print_top_combined(rankings, top_n=10)
    print_best_recommendations(rankings, top_n=5)

    print(f"\nFull result saved to: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()