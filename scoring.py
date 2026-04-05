from __future__ import annotations

import math
from typing import Any, Dict, List, Optional


# ============================================================
# Tunable constants
# ============================================================

# Expression scoring
EXPR_MAX_SCORE = 35.0
NEAR_RATIO_MAX_SCORE = 45.0
FAR_RATIO_MAX_SCORE = 15.0
DOMINANT_TISSUE_BONUS = 5.0

# Ratio / expression transforms
EXPR_REF = 50.0                  # reference TPM-like scale for saturation
RATIO_PSEUDOCOUNT = 0.1          # protects against divide-by-zero / sparse RNA-seq
RATIO_SATURATION = 100.0         # ratio at which ratio-score saturates
NEAR_SCALE_REF_EXPR = 10.0       # controls how much target_expr gates near-score
FAR_SCALE_REF_EXPR = 10.0        # controls how much target_expr gates far-score

# Detection thresholds
VERY_LOW_EXPR_THRESHOLD = 1.0
LOW_EXPR_THRESHOLD = 5.0
DETECTION_FLOOR = 0.5            # below this, keep but mark low-confidence
VERY_WEAK_FLOOR = 0.2            # even weaker

# Weak-ratio warnings
WEAK_RATIO_THRESHOLD = 1.5

# Combined ranking mild adjustments
THERMO_ADJUSTMENTS = {
    "excellent": 2.0,
    "good": 1.0,
    "usable": 0.0,
    "warning": -2.0,
}
SPECIFICITY_ADJUSTMENTS = {
    "clean": 2.0,
    "far_off_target_warning": -1.0,
    "near_off_target_warning": -4.0,
    "failed_target": -20.0,
}


# ============================================================
# Public API
# ============================================================

def build_rankings(
    d_results: List[Dict[str, Any]],
    e_evals: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Build three rankings:
    1. thermo_ranking     -> pair-level
    2. expression_ranking -> pair+source-level
    3. combined_ranking   -> pair+source-level
    """
    e_map = {item.get("pair_index"): item for item in e_evals}

    pair_profiles: List[Dict[str, Any]] = []
    pair_source_profiles: List[Dict[str, Any]] = []

    for d_result in d_results:
        pair_index = d_result.get("pair_index")
        e_eval = e_map.get(pair_index, {})

        pair_profile = build_pair_profile(d_result)
        pair_profiles.append(pair_profile)

        source_profiles = build_pair_source_profiles(d_result, e_eval)
        pair_source_profiles.extend(source_profiles)

    thermo_ranking = rank_thermo(pair_profiles)
    expression_ranking = rank_expression(pair_source_profiles)
    combined_ranking = rank_combined(pair_profiles, pair_source_profiles)

    return {
        "thermo_ranking": thermo_ranking,
        "expression_ranking": expression_ranking,
        "combined_ranking": combined_ranking,
    }


def build_pair_profile(d_result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Build pair-level profile:
    - thermo label
    - specificity label
    - warnings
    - thermo score (for thermo ranking only)
    """
    pair_index = d_result.get("pair_index")

    forward = d_result.get("forward")
    reverse = d_result.get("reverse")

    tm_diff = _to_float(d_result.get("tm_diff"), default=0.0)
    heterodimer_dG = _to_float(d_result.get("heterodimer_dG"), default=0.0)
    product_size = d_result.get("product_size")

    target_amplicons = d_result.get("target_amplicons", []) or []
    off_targets = d_result.get("off_targets", {}) or {}
    near_hits = off_targets.get("near", []) or []
    far_hits = off_targets.get("far", []) or []

    thermo_warnings: List[str] = []
    specificity_warnings: List[str] = []

    thermo_score = _compute_thermo_score(
        tm_diff=tm_diff,
        heterodimer_dG=heterodimer_dG,
        product_size=product_size,
        warnings=thermo_warnings,
    )
    thermo_label = _classify_thermo_label(
        tm_diff=tm_diff,
        heterodimer_dG=heterodimer_dG,
        product_size=product_size,
    )

    specificity_label = _classify_specificity_label(
        target_amplicon_count=len(target_amplicons),
        near_off_target_count=len(near_hits),
        far_off_target_count=len(far_hits),
        warnings=specificity_warnings,
    )

    return {
        "pair_index": pair_index,
        "forward": _extract_primer_sequence(forward),
        "reverse": _extract_primer_sequence(reverse),
        "thermo": {
            "tm_diff": tm_diff,
            "heterodimer_dG": heterodimer_dG,
            "product_size": product_size,
            "thermo_score": round(thermo_score, 3),
            "label": thermo_label,
            "warnings": thermo_warnings,
        },
        "specificity": {
            "target_amplicon_count": len(target_amplicons),
            "near_off_target_count": len(near_hits),
            "far_off_target_count": len(far_hits),
            "label": specificity_label,
            "warnings": specificity_warnings,
        },
    }


def build_pair_source_profiles(
    d_result: Dict[str, Any],
    e_eval: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Build pair+source(tissue)-level profiles.
    Recommendation unit:
    primer pair + PCR source/tissue
    """
    pair_index = d_result.get("pair_index")

    target_expr_by_tissue = e_eval.get("target_expr_by_tissue", {}) or {}
    near_expr_by_tissue = e_eval.get("near_expr_by_tissue", {}) or {}
    far_expr_by_tissue = e_eval.get("far_expr_by_tissue", {}) or {}

    target_vs_near_ratio = e_eval.get("target_vs_near_ratio", {}) or {}
    target_vs_far_ratio = e_eval.get("target_vs_far_ratio", {}) or {}
    dominant_target_tissues = e_eval.get("dominant_target_tissues", []) or []

    tissues = _collect_tissues(
        target_expr_by_tissue,
        near_expr_by_tissue,
        far_expr_by_tissue,
        target_vs_near_ratio,
        target_vs_far_ratio,
    )

    profiles: List[Dict[str, Any]] = []

    for tissue in tissues:
        target_expr = _safe_lookup(target_expr_by_tissue, tissue)
        near_expr = _safe_lookup(near_expr_by_tissue, tissue)
        far_expr = _safe_lookup(far_expr_by_tissue, tissue)

        near_ratio_raw = _safe_lookup_or_none(target_vs_near_ratio, tissue)
        far_ratio_raw = _safe_lookup_or_none(target_vs_far_ratio, tissue)

        near_ratio = _normalize_ratio(
            ratio_value=near_ratio_raw,
            target_expr=target_expr,
            competitor_expr=near_expr,
        )
        far_ratio = _normalize_ratio(
            ratio_value=far_ratio_raw,
            target_expr=target_expr,
            competitor_expr=far_expr,
        )

        score_detail = _compute_expression_score_detail(
            target_expr=target_expr,
            near_expr=near_expr,
            far_expr=far_expr,
            target_vs_near_ratio=near_ratio,
            target_vs_far_ratio=far_ratio,
            is_dominant_tissue=(tissue in dominant_target_tissues),
        )
        expression_score = score_detail["expression_score"]

        recommendation_label = _classify_expression_label(
            expression_score=expression_score,
            target_expr=target_expr,
            near_expr=near_expr,
            far_expr=far_expr,
        )

        source_warnings = _build_source_warnings(
            tissue=tissue,
            target_expr=target_expr,
            near_expr=near_expr,
            far_expr=far_expr,
            target_vs_near_ratio=near_ratio,
            target_vs_far_ratio=far_ratio,
            is_dominant_tissue=(tissue in dominant_target_tissues),
        )

        source_tags = _build_source_tags(
            target_expr=target_expr,
            near_expr=near_expr,
            far_expr=far_expr,
            expression_score=expression_score,
            is_dominant_tissue=(tissue in dominant_target_tissues),
        )

        profiles.append({
            "pair_index": pair_index,
            "source": tissue,
            "target_expr": round(target_expr, 6),
            "near_expr": round(near_expr, 6),
            "far_expr": round(far_expr, 6),
            "target_vs_near_ratio": round(near_ratio, 6),
            "target_vs_far_ratio": round(far_ratio, 6),
            "expression_score": round(expression_score, 3),
            "recommendation_label": recommendation_label,
            "is_dominant_target_tissue": tissue in dominant_target_tissues,
            "source_tags": source_tags,
            "source_warnings": source_warnings,
            "score_detail": {
                "expr_component": round(score_detail["expr_component"], 3),
                "near_ratio_component": round(score_detail["near_ratio_component"], 3),
                "far_ratio_component": round(score_detail["far_ratio_component"], 3),
                "dominant_bonus": round(score_detail["dominant_bonus"], 3),
                "low_expression_penalty": round(score_detail["low_expression_penalty"], 3),
                "near_competition_penalty": round(score_detail["near_competition_penalty"], 3),
                "far_competition_penalty": round(score_detail["far_competition_penalty"], 3),
            },
        })

    return profiles


def rank_thermo(pair_profiles: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []

    for pair in pair_profiles:
        thermo = pair.get("thermo", {})
        specificity = pair.get("specificity", {})

        rows.append({
            "pair_index": pair.get("pair_index"),
            "forward": pair.get("forward"),
            "reverse": pair.get("reverse"),
            "thermo_score": thermo.get("thermo_score", 0.0),
            "thermo_label": thermo.get("label"),
            "specificity_label": specificity.get("label"),
            "warnings": [
                *thermo.get("warnings", []),
                *specificity.get("warnings", []),
            ],
        })

    rows.sort(key=lambda x: x["thermo_score"], reverse=True)
    return rows


def rank_expression(pair_source_profiles: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = [dict(item) for item in pair_source_profiles]
    rows.sort(key=lambda x: x["expression_score"], reverse=True)
    return rows


def rank_combined(
    pair_profiles: List[Dict[str, Any]],
    pair_source_profiles: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    pair_map = {item.get("pair_index"): item for item in pair_profiles}
    rows: List[Dict[str, Any]] = []

    for item in pair_source_profiles:
        pair_index = item.get("pair_index")
        pair_profile = pair_map.get(pair_index, {})

        thermo = pair_profile.get("thermo", {})
        specificity = pair_profile.get("specificity", {})

        expression_score = _to_float(item.get("expression_score"), default=0.0)
        thermo_score = _to_float(thermo.get("thermo_score"), default=0.0)

        thermo_label = thermo.get("label", "unknown")
        specificity_label = specificity.get("label", "unknown")

        combined_score = expression_score
        combined_score += _thermo_adjustment(thermo_label, thermo_score)
        combined_score += _specificity_adjustment(specificity_label)
        combined_score = max(0.0, min(100.0, combined_score))

        final_warnings = [
            *item.get("source_warnings", []),
            *thermo.get("warnings", []),
            *specificity.get("warnings", []),
        ]

        final_tags = [
            *item.get("source_tags", []),
            thermo_label,
            specificity_label,
        ]

        rows.append({
            "pair_index": pair_index,
            "source": item.get("source"),
            "combined_score": round(combined_score, 3),
            "expression_score": round(expression_score, 3),
            "thermo_score": round(thermo_score, 3),
            "thermo_label": thermo_label,
            "specificity_label": specificity_label,
            "target_expr": item.get("target_expr"),
            "near_expr": item.get("near_expr"),
            "far_expr": item.get("far_expr"),
            "target_vs_near_ratio": item.get("target_vs_near_ratio"),
            "target_vs_far_ratio": item.get("target_vs_far_ratio"),
            "final_recommendation": _classify_combined_label(combined_score),
            "source_tags": item.get("source_tags", []),
            "source_warnings": item.get("source_warnings", []),
            "all_tags": _dedup_preserve_order(final_tags),
            "warnings": _dedup_preserve_order(final_warnings),
            "score_detail": item.get("score_detail", {}),
        })

    rows.sort(key=lambda x: x["combined_score"], reverse=True)
    return rows


def select_top_recommendations(
    rankings: Dict[str, List[Dict[str, Any]]],
    *,
    ranking_key: str = "combined_ranking",
    top_n: int = 10,
) -> List[Dict[str, Any]]:
    rows = rankings.get(ranking_key, []) or []
    return rows[:top_n]


# ============================================================
# Core scoring logic
# ============================================================

def _compute_thermo_score(
    *,
    tm_diff: float,
    heterodimer_dG: float,
    product_size: Any,
    warnings: List[str],
) -> float:
    score = 100.0

    if tm_diff > 1.0:
        score -= min((tm_diff - 1.0) * 10.0, 25.0)
        warnings.append(f"Tm difference slightly high: {tm_diff:.2f}")
    if tm_diff > 2.0:
        warnings.append("Tm difference exceeds ideal range.")

    if heterodimer_dG <= -9.0:
        score -= 25.0
        warnings.append("Strong heterodimer risk.")
    elif heterodimer_dG <= -6.0:
        score -= 12.0
        warnings.append("Moderate heterodimer risk.")

    if isinstance(product_size, (int, float)):
        if product_size < 80 or product_size > 1500:
            score -= 8.0
            warnings.append("Product size outside preferred range (80-1500 bp).")

    return max(0.0, min(100.0, score))


def _compute_expression_score_detail(
    *,
    target_expr: float,
    near_expr: float,
    far_expr: float,
    target_vs_near_ratio: float,
    target_vs_far_ratio: float,
    is_dominant_tissue: bool,
) -> Dict[str, float]:
    """
    New scoring model:
    - expression uses log-scale
    - ratio uses pseudo-count + log-scale
    - near/far ratio contribution is gated by target_expr
    - low-expression tissues are kept, but penalized / tagged
    """
    expr_component = _bounded_expr_score_log(target_expr, max_score=EXPR_MAX_SCORE)

    near_ratio_component = _bounded_ratio_score_log(
        target_vs_near_ratio,
        max_score=NEAR_RATIO_MAX_SCORE,
        saturation_ratio=RATIO_SATURATION,
    )
    far_ratio_component = _bounded_ratio_score_log(
        target_vs_far_ratio,
        max_score=FAR_RATIO_MAX_SCORE,
        saturation_ratio=RATIO_SATURATION,
    )

    # Gate near/far contributions by target expression strength.
    near_scale = _target_expr_scale(target_expr, ref_expr=NEAR_SCALE_REF_EXPR)
    far_scale = _target_expr_scale(target_expr, ref_expr=FAR_SCALE_REF_EXPR)

    near_ratio_component *= near_scale
    far_ratio_component *= far_scale

    dominant_bonus = DOMINANT_TISSUE_BONUS if is_dominant_tissue else 0.0

    low_expression_penalty = 0.0
    if 0 < target_expr < DETECTION_FLOOR:
        low_expression_penalty = 12.0
    elif DETECTION_FLOOR <= target_expr < VERY_LOW_EXPR_THRESHOLD:
        low_expression_penalty = 6.0
    elif VERY_LOW_EXPR_THRESHOLD <= target_expr < LOW_EXPR_THRESHOLD:
        low_expression_penalty = 2.0

    near_competition_penalty = 0.0
    far_competition_penalty = 0.0

    if near_expr > 0 and target_expr <= near_expr:
        near_competition_penalty = 10.0

    if far_expr > 0 and target_expr <= far_expr:
        far_competition_penalty = 5.0

    expression_score = (
        expr_component
        + near_ratio_component
        + far_ratio_component
        + dominant_bonus
        - low_expression_penalty
        - near_competition_penalty
        - far_competition_penalty
    )
    expression_score = max(0.0, min(100.0, expression_score))

    return {
        "expression_score": expression_score,
        "expr_component": expr_component,
        "near_ratio_component": near_ratio_component,
        "far_ratio_component": far_ratio_component,
        "dominant_bonus": dominant_bonus,
        "low_expression_penalty": low_expression_penalty,
        "near_competition_penalty": near_competition_penalty,
        "far_competition_penalty": far_competition_penalty,
    }


# ============================================================
# Labels / tags / warnings
# ============================================================

def _classify_thermo_label(
    *,
    tm_diff: float,
    heterodimer_dG: float,
    product_size: Any,
) -> str:
    if tm_diff <= 1.0 and heterodimer_dG > -6.0 and _is_preferred_product_size(product_size):
        return "excellent"
    if tm_diff <= 2.0 and heterodimer_dG > -9.0:
        return "good"
    if tm_diff <= 3.0:
        return "usable"
    return "warning"


def _classify_specificity_label(
    *,
    target_amplicon_count: int,
    near_off_target_count: int,
    far_off_target_count: int,
    warnings: List[str],
) -> str:
    if target_amplicon_count <= 0:
        warnings.append("No target amplicon detected.")
        return "failed_target"

    if near_off_target_count == 0 and far_off_target_count == 0:
        return "clean"

    if near_off_target_count > 0:
        warnings.append(f"{near_off_target_count} near off-target amplicon(s).")
        if near_off_target_count >= 3:
            warnings.append("Multiple near off-targets may affect interpretation.")
        return "near_off_target_warning"

    warnings.append(f"{far_off_target_count} far off-target amplicon(s).")
    return "far_off_target_warning"


def _classify_expression_label(
    *,
    expression_score: float,
    target_expr: float,
    near_expr: float,
    far_expr: float,
) -> str:
    if target_expr <= 0:
        return "not_expressed"

    if target_expr < VERY_WEAK_FLOOR:
        return "very_low_confidence"

    if expression_score >= 85:
        return "highly_recommended"
    if expression_score >= 70:
        return "recommended"
    if expression_score >= 50:
        return "usable"

    if near_expr > target_expr:
        return "near_competition_risk"

    if far_expr > target_expr:
        return "far_competition_risk"

    return "weak"


def _classify_combined_label(score: float) -> str:
    if score >= 85:
        return "excellent"
    if score >= 70:
        return "good"
    if score >= 50:
        return "usable"
    return "warning"


def _build_source_tags(
    *,
    target_expr: float,
    near_expr: float,
    far_expr: float,
    expression_score: float,
    is_dominant_tissue: bool,
) -> List[str]:
    tags: List[str] = []

    if is_dominant_tissue:
        tags.append("dominant_target_tissue")

    if target_expr <= 0:
        tags.append("not_expressed")
    elif target_expr < VERY_WEAK_FLOOR:
        tags.append("very_low_confidence")
    elif target_expr < VERY_LOW_EXPR_THRESHOLD:
        tags.append("very_low_expression")
    elif target_expr < LOW_EXPR_THRESHOLD:
        tags.append("low_expression")
    else:
        tags.append("expressed")

    if near_expr > 0 and target_expr <= near_expr:
        tags.append("near_competition")
    if far_expr > 0 and target_expr <= far_expr:
        tags.append("far_competition")

    if expression_score >= 85:
        tags.append("expression_top")
    elif expression_score >= 70:
        tags.append("expression_good")

    return tags


def _build_source_warnings(
    *,
    tissue: str,
    target_expr: float,
    near_expr: float,
    far_expr: float,
    target_vs_near_ratio: float,
    target_vs_far_ratio: float,
    is_dominant_tissue: bool,
) -> List[str]:
    warnings: List[str] = []

    if target_expr <= 0:
        warnings.append(f"{tissue}: target not detectably expressed.")
    elif target_expr < VERY_WEAK_FLOOR:
        warnings.append(
            f"{tissue}: target expression is extremely low; keep in list but treat as low-confidence for PCR."
        )
    elif target_expr < DETECTION_FLOOR:
        warnings.append(
            f"{tissue}: target expression is below detection-favorable range; interpret cautiously."
        )
    elif target_expr < VERY_LOW_EXPR_THRESHOLD:
        warnings.append(f"{tissue}: target expression is very low; keep in list but interpret cautiously.")
    elif target_expr < LOW_EXPR_THRESHOLD:
        warnings.append(f"{tissue}: target expression is low.")

    if near_expr > 0 and target_expr <= near_expr:
        warnings.append(f"{tissue}: near off-target expression competes strongly with target.")

    if far_expr > 0 and target_expr <= far_expr:
        warnings.append(f"{tissue}: far off-target expression may interfere with interpretation.")

    if near_expr > 0 and target_vs_near_ratio < WEAK_RATIO_THRESHOLD:
        warnings.append(f"{tissue}: target/near ratio is weak.")

    if far_expr > 0 and target_vs_far_ratio < WEAK_RATIO_THRESHOLD:
        warnings.append(f"{tissue}: target/far ratio is weak.")

    if is_dominant_tissue and target_expr < VERY_LOW_EXPR_THRESHOLD:
        warnings.append(f"{tissue}: marked dominant by E, but absolute target expression is still very low.")

    return warnings


# ============================================================
# Combined ranking adjustments
# ============================================================

def _thermo_adjustment(thermo_label: str, thermo_score: float) -> float:
    return THERMO_ADJUSTMENTS.get(thermo_label, 0.0)


def _specificity_adjustment(specificity_label: str) -> float:
    return SPECIFICITY_ADJUSTMENTS.get(specificity_label, 0.0)


# ============================================================
# Helper functions
# ============================================================

def _extract_primer_sequence(primer_obj: Any) -> Optional[str]:
    if isinstance(primer_obj, dict):
        return primer_obj.get("sequence")
    if isinstance(primer_obj, str):
        return primer_obj
    return None


def _collect_tissues(*dicts: Dict[str, Any]) -> List[str]:
    tissues = set()
    for d in dicts:
        if isinstance(d, dict):
            tissues.update(d.keys())
    return sorted(tissues)


def _safe_lookup(d: Dict[str, Any], key: str, default: float = 0.0) -> float:
    value = d.get(key, default)
    return _to_float(value, default=default)


def _safe_lookup_or_none(d: Dict[str, Any], key: str) -> Optional[float]:
    if key not in d:
        return None
    value = d.get(key)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _normalize_ratio(
    *,
    ratio_value: Optional[float],
    target_expr: float,
    competitor_expr: float,
) -> float:
    """
    New ratio behavior:
    - If ratio already exists, use it
    - Else compute from pseudo-count-smoothed expression
    """
    if ratio_value is not None:
        return max(0.0, float(ratio_value))

    # Pseudo-count smoothing prevents infinite / zero-collapse behavior.
    return (target_expr + RATIO_PSEUDOCOUNT) / (competitor_expr + RATIO_PSEUDOCOUNT)


def _bounded_expr_score_log(expr: float, *, max_score: float) -> float:
    """
    Log-scaled expression score.
    Saturates around EXPR_REF.
    """
    if expr <= 0:
        return 0.0

    numerator = math.log10(expr + 1.0)
    denominator = math.log10(EXPR_REF + 1.0)
    if denominator <= 0:
        return 0.0

    frac = min(1.0, numerator / denominator)
    return max_score * frac


def _bounded_ratio_score_log(
    ratio: float,
    *,
    max_score: float,
    saturation_ratio: float,
) -> float:
    """
    Log-scaled ratio score.
    ratio <= 1 gives weak signal
    ratio around saturation_ratio gives full score.
    """
    if ratio <= 0:
        return 0.0

    numerator = math.log10(ratio + 1.0)
    denominator = math.log10(saturation_ratio + 1.0)
    if denominator <= 0:
        return 0.0

    frac = min(1.0, numerator / denominator)
    return max_score * frac


def _target_expr_scale(expr: float, *, ref_expr: float) -> float:
    """
    Gates ratio contribution by target expression strength.
    When target expression is extremely low, even a 'clean' ratio should not dominate.
    """
    if expr <= 0:
        return 0.0

    numerator = math.log10(expr + 1.0)
    denominator = math.log10(ref_expr + 1.0)
    if denominator <= 0:
        return 0.0

    return min(1.0, numerator / denominator)


def _is_preferred_product_size(product_size: Any) -> bool:
    if not isinstance(product_size, (int, float)):
        return False
    return 80 <= product_size <= 1500


def _dedup_preserve_order(items: List[str]) -> List[str]:
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out