from __future__ import annotations

import re
from typing import Any, Callable, Dict, List, Optional, Tuple

from scan_primers import scan_primers
from filter_primers import filter_primers
from pair_builder import build_pairs
from blast_client_d import process_c_module_output_d

from E.expression.client import ExpressionClient
from E.evaluate import evaluate_blast_results
from E.resolver.resolver import TranscriptResolver

from scoring import build_rankings


def _clean_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    seq = re.sub(r"\s+", "", sequence).upper()
    if not seq:
        raise ValueError("sequence is empty after cleaning")

    bad = re.findall(r"[^ACGT]", seq)
    if bad:
        raise ValueError(f"sequence contains invalid characters: {sorted(set(bad))}")

    return seq


def _build_resolver_fn(
    resolver_obj: TranscriptResolver,
) -> Callable[[Dict[str, Any]], Optional[str]]:
    def _looks_like_refseq(x: str) -> bool:
        if not x:
            return False
        return x.startswith(("XM_", "NM_", "XR_", "NR_"))

    def _looks_like_glyma(x: str) -> bool:
        if not x:
            return False
        return x.startswith("Glyma.")

    def _resolver_fn(amplicon: Dict[str, Any]) -> Optional[str]:
        glyma_id = (amplicon.get("glyma_id") or "").strip()
        if glyma_id:
            return glyma_id

        gene_id = (amplicon.get("gene_id") or "").strip()
        if _looks_like_glyma(gene_id):
            return gene_id

        subject_id = (amplicon.get("subject_id") or "").strip()
        if _looks_like_refseq(subject_id):
            rr = resolver_obj.resolve_one(subject_id)
            gid = (rr.get("glyma_id") or "").strip()
            if gid:
                return gid

        if _looks_like_refseq(gene_id):
            rr = resolver_obj.resolve_one(gene_id)
            gid = (rr.get("glyma_id") or "").strip()
            if gid:
                return gid

        return None

    return _resolver_fn


def _collect_real_expression_sources(
    evaluated_blast_results: List[Dict[str, Any]],
) -> List[str]:
    source_set = set()

    for pair in evaluated_blast_results:
        ee = pair.get("e_eval") or {}
        expr_src = ee.get("expression_source_detail") or {}

        for bucket in ("target", "near", "far"):
            bucket_detail = expr_src.get(bucket) or {}
            for _gid, detail in bucket_detail.items():
                if not isinstance(detail, dict):
                    continue

                src = detail.get("source")
                if src:
                    source_set.add(str(src))

                by_source = detail.get("by_source") or {}
                if isinstance(by_source, dict):
                    for bs in by_source.keys():
                        source_set.add(str(bs))

    return sorted(source_set)


def _build_scoring_inputs(
    blast_results: List[Dict[str, Any]],
    evaluated_blast_results: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    将 D / E 的输出整理成 scoring.py 需要的输入格式。

    返回:
        d_results_for_scoring, e_evals_for_scoring
    """
    d_results_for_scoring: List[Dict[str, Any]] = []
    e_evals_for_scoring: List[Dict[str, Any]] = []

    eval_map: Dict[Any, Dict[str, Any]] = {
        item.get("pair_index"): item for item in evaluated_blast_results
    }

    for idx, d_item in enumerate(blast_results):
        pair_index = d_item.get("pair_index", idx)

        d_result = dict(d_item)
        d_result["pair_index"] = pair_index
        d_results_for_scoring.append(d_result)

        evaluated = eval_map.get(pair_index, {}) or {}
        raw_e_eval = evaluated.get("e_eval") or {}

        e_eval = {
            "pair_index": pair_index,
            "target_expr_by_tissue": raw_e_eval.get("target_expr_by_tissue", {}) or {},
            "near_expr_by_tissue": raw_e_eval.get("near_expr_by_tissue", {}) or {},
            "far_expr_by_tissue": raw_e_eval.get("far_expr_by_tissue", {}) or {},
            "target_vs_near_ratio": raw_e_eval.get("target_vs_near_ratio", {}) or {},
            "target_vs_far_ratio": raw_e_eval.get("target_vs_far_ratio", {}) or {},
            "dominant_target_tissues": raw_e_eval.get("dominant_target_tissues", []) or [],
        }
        e_evals_for_scoring.append(e_eval)

    return d_results_for_scoring, e_evals_for_scoring


def pipeline_ABCDES(
    sequence: str,
    # A
    size_range: Tuple[int, int] = (18, 25),
    # B
    min_tm: float = 50.0,
    max_tm: float = 65.0,
    max_hairpin_dG: float = -3.0,
    max_dimer_dG: float = -9.0,
    min_gc: float = 0.0,
    max_gc: float = 1.0,
    min_len: int = 18,
    max_len: int = 25,
    # C
    target_length: Optional[int] = None,
    max_tm_diff: float = 2.0,
    max_heterodimer_dG: float = 0.0,
    # D
    organism: str = "Glycine max",
    db: str = "refseq_rna",
    max_workers: int = 3,
    sleep_between_requests: float = 0.5,
    require_strict_orientation: bool = True,
    min_amplicon_len: int = 50,
    max_amplicon_len: int = 3000,
    near_bp: int = 30,
    target_subject_id: Optional[str] = None,
    dummy_blast: bool = False,
    min_identity_min: Optional[float] = None,
    # E
    enable_e: bool = True,
    expression_cache_db_path: str = "expression_cache.sqlite",
    expression_ttl_seconds: int = 7 * 24 * 3600,
    allowed_canon_conditions: Optional[set[str]] = None,
    conekt_prefix: str = "glyma.Wm82.gnm2.ann1.",
    conekt_sleep_s: float = 0.2,
    conekt_timeout: int = 25,
    e_use_cache: bool = True,
    preferred_expression_source: str = "integrated",
    enable_soyomics_fallback: bool = True,
    soyomics_timeout: int = 25,
    source_weights: Optional[Dict[str, float]] = None,
    # resolver
    resolver_cache_db_path: str = "resolver_cache.sqlite",
    enable_ensembl_fallback: bool = True,
    # optional external resolver
    e_resolver: Optional[Callable[[Dict[str, Any]], Optional[str]]] = None,
    # S
    enable_scoring: bool = True,
) -> Dict[str, Any]:
    seq = _clean_sequence(sequence)
    if target_length is None:
        target_length = len(seq)

    # A
    scan_success, scan_result, scan_status = scan_primers(seq, size_range=size_range)
    if not scan_success:
        return {
            "scan": scan_result,
            "filter": {"forward": [], "reverse": []},
            "pairs": [],
            "blast_results": [],
            "evaluated_blast_results": [],
            "rankings": {
                "thermo_ranking": [],
                "expression_ranking": [],
                "combined_ranking": [],
            },
            "meta": {
                "pipeline": "ABCDES",
                "stage_failed": "A",
                "error": scan_status.get("error", "scan_primers failed"),
            },
        }

    forward_primers = scan_result.get("forward", [])
    reverse_primers = scan_result.get("reverse", [])

    # B
    filtered_forward = filter_primers(
        forward_primers,
        min_tm=min_tm,
        max_tm=max_tm,
        max_hairpin_dG=max_hairpin_dG,
        max_dimer_dG=max_dimer_dG,
        min_gc=min_gc,
        max_gc=max_gc,
        min_len=min_len,
        max_len=max_len,
    )

    filtered_reverse = filter_primers(
        reverse_primers,
        min_tm=min_tm,
        max_tm=max_tm,
        max_hairpin_dG=max_hairpin_dG,
        max_dimer_dG=max_dimer_dG,
        min_gc=min_gc,
        max_gc=max_gc,
        min_len=min_len,
        max_len=max_len,
    )

    # C
    pairs = build_pairs(
        filtered_forward,
        filtered_reverse,
        target_length=target_length,
        max_tm_diff=max_tm_diff,
        max_heterodimer_dG=max_heterodimer_dG,
    )

    # D
    blast_results = process_c_module_output_d(
        c_module_output=pairs,
        organism=organism,
        db=db,
        max_workers=max_workers,
        target_length=target_length,
        near_bp=near_bp,
        target_subject_id=target_subject_id,
        dummy=dummy_blast,
        min_amplicon_len=min_amplicon_len,
        max_amplicon_len=max_amplicon_len,
        require_strict_orientation=require_strict_orientation,
        sleep_between_requests=sleep_between_requests,
        min_identity_min=min_identity_min,
    )

    # E
    evaluated_blast_results: List[Dict[str, Any]] = []
    e_status: Dict[str, Any] = {
        "enabled": enable_e,
        "ok": False,
        "error": None,
        "resolver_enabled": False,
        "resolved_target_pairs": 0,
        "pairs_with_nonempty_target_expr": 0,
        "expression_sources_used": [],
    }

    if enable_e:
        try:
            expression_client = ExpressionClient(
                cache_db_path=expression_cache_db_path,
                ttl_seconds=expression_ttl_seconds,
                allowed_canon_conditions=allowed_canon_conditions,
                conekt_prefix=conekt_prefix,
                conekt_sleep_s=conekt_sleep_s,
                conekt_timeout=conekt_timeout,
                preferred_source=preferred_expression_source,
                enable_soyomics_fallback=enable_soyomics_fallback,
                soyomics_timeout=soyomics_timeout,
                source_weights=source_weights,
            )

            if e_resolver is not None:
                resolver_fn = e_resolver
                e_status["resolver_enabled"] = True
            else:
                transcript_resolver = TranscriptResolver(
                    cache_path=resolver_cache_db_path,
                    enable_ensembl_fallback=enable_ensembl_fallback,
                )
                resolver_fn = _build_resolver_fn(transcript_resolver)
                e_status["resolver_enabled"] = True

            evaluated_blast_results = evaluate_blast_results(
                blast_results,
                expression_client=expression_client,
                resolver=resolver_fn,
                use_cache=e_use_cache,
                attach_key="e_eval",
            )

            resolved_target_pairs = 0
            pairs_with_nonempty_target_expr = 0

            for pair in evaluated_blast_results:
                ee = pair.get("e_eval") or {}
                if ee.get("target_genes"):
                    resolved_target_pairs += 1
                if ee.get("target_expr_by_tissue"):
                    pairs_with_nonempty_target_expr += 1

            e_status["resolved_target_pairs"] = resolved_target_pairs
            e_status["pairs_with_nonempty_target_expr"] = pairs_with_nonempty_target_expr
            e_status["expression_sources_used"] = _collect_real_expression_sources(
                evaluated_blast_results
            )
            e_status["ok"] = True

        except Exception as e:
            evaluated_blast_results = [dict(x) for x in blast_results]
            e_status["ok"] = False
            e_status["error"] = f"{type(e).__name__}: {e}"

    else:
        evaluated_blast_results = [dict(x) for x in blast_results]

    # S
    rankings: Dict[str, Any] = {
        "thermo_ranking": [],
        "expression_ranking": [],
        "combined_ranking": [],
    }
    scoring_status: Dict[str, Any] = {
        "enabled": enable_scoring,
        "ok": False,
        "error": None,
        "n_pairs_in_d": 0,
        "n_pairs_in_e": 0,
    }

    if enable_scoring:
        try:
            d_results_for_scoring, e_evals_for_scoring = _build_scoring_inputs(
                blast_results=blast_results,
                evaluated_blast_results=evaluated_blast_results,
            )

            rankings = build_rankings(
                d_results=d_results_for_scoring,
                e_evals=e_evals_for_scoring,
            )

            scoring_status["ok"] = True
            scoring_status["n_pairs_in_d"] = len(d_results_for_scoring)
            scoring_status["n_pairs_in_e"] = len(e_evals_for_scoring)

        except Exception as e:
            scoring_status["ok"] = False
            scoring_status["error"] = f"{type(e).__name__}: {e}"

    meta = {
        "pipeline": "ABCDES",
        "target_length": target_length,
        "organism": organism,
        "db": db,
        "d_params": {
            "max_workers": max_workers,
            "sleep_between_requests": sleep_between_requests,
            "require_strict_orientation": require_strict_orientation,
            "min_amplicon_len": min_amplicon_len,
            "max_amplicon_len": max_amplicon_len,
            "near_bp": near_bp,
            "target_subject_id": target_subject_id,
            "dummy_blast": dummy_blast,
            "min_identity_min": min_identity_min,
        },
        "e_params": {
            "enable_e": enable_e,
            "expression_cache_db_path": expression_cache_db_path,
            "expression_ttl_seconds": expression_ttl_seconds,
            "allowed_canon_conditions": sorted(allowed_canon_conditions)
            if allowed_canon_conditions
            else None,
            "conekt_prefix": conekt_prefix,
            "conekt_sleep_s": conekt_sleep_s,
            "conekt_timeout": conekt_timeout,
            "e_use_cache": e_use_cache,
            "preferred_expression_source": preferred_expression_source,
            "enable_soyomics_fallback": enable_soyomics_fallback,
            "soyomics_timeout": soyomics_timeout,
            "source_weights": source_weights or {},
        },
        "resolver_params": {
            "resolver_cache_db_path": resolver_cache_db_path,
            "enable_ensembl_fallback": enable_ensembl_fallback,
            "custom_resolver_provided": e_resolver is not None,
        },
        "e_status": e_status,
        "scoring_status": scoring_status,
        "counts": {
            "forward_primers": len(forward_primers),
            "reverse_primers": len(reverse_primers),
            "filtered_forward": len(filtered_forward),
            "filtered_reverse": len(filtered_reverse),
            "pairs": len(pairs),
            "blast_results": len(blast_results),
            "evaluated_blast_results": len(evaluated_blast_results),
            "thermo_ranking": len(rankings.get("thermo_ranking", []) or []),
            "expression_ranking": len(rankings.get("expression_ranking", []) or []),
            "combined_ranking": len(rankings.get("combined_ranking", []) or []),
        },
    }

    return {
        "scan": scan_result,
        "filter": {
            "forward": filtered_forward,
            "reverse": filtered_reverse,
        },
        "pairs": pairs,
        "blast_results": blast_results,
        "evaluated_blast_results": evaluated_blast_results,
        "rankings": rankings,
        "meta": meta,
    }