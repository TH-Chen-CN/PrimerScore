from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional


EPS = 1e-9


def _top_tissues(expr: Dict[str, float], n: int = 3) -> List[str]:
    """
    只返回 > 0 的 top tissues。
    避免以前那种全是 0 还把 flower/leaf 排进来的情况。
    """
    items = [(k, float(v)) for k, v in expr.items() if float(v) > 0]
    items.sort(key=lambda x: x[1], reverse=True)
    return [k for k, _ in items[:n]]


def _merge_expr_for_gene_ids(
    gene_ids: List[str],
    expression_client: Any,
    use_cache: bool = True,
) -> Dict[str, Any]:
    """
    聚合多个 gene 的表达。
    当前策略：
    - 每个 tissue 取 max
    - raw 细节全部保留
    """
    expr_by_tissue: Dict[str, float] = {}
    raw_detail: Dict[str, Dict[str, float]] = {}
    raw_to_canonical_detail: Dict[str, Dict[str, Optional[str]]] = {}
    missing_expression_genes: List[str] = []
    source_details: Dict[str, Any] = {}

    for gid in gene_ids:
        ge = expression_client.fetch_gene(gid, use_cache=use_cache)

        source_details[gid] = {
            "source": ge.source,
            "notes": ge.notes,
            "source_url": ge.source_url,
            "by_source": ge.by_source,
        }

        if not ge.ok:
            missing_expression_genes.append(gid)
            continue

        for tissue, value in ge.canonical_expr.items():
            expr_by_tissue[tissue] = max(float(value), float(expr_by_tissue.get(tissue, 0.0)))

        if ge.raw_expr:
            raw_detail[gid] = dict(ge.raw_expr)
        if ge.raw_to_canonical:
            raw_to_canonical_detail[gid] = dict(ge.raw_to_canonical)

    return {
        "expr_by_tissue": expr_by_tissue,
        "raw_detail": raw_detail,
        "raw_to_canonical_detail": raw_to_canonical_detail,
        "missing_expression_genes": missing_expression_genes,
        "source_details": source_details,
    }


def _resolve_amplicon_gene_ids(
    amplicons: List[Dict[str, Any]],
    resolver: Optional[Callable[[Dict[str, Any]], Optional[str]]] = None,
) -> Dict[str, Any]:
    gene_ids: List[str] = []
    missing: List[Dict[str, Any]] = []

    for amp in amplicons:
        gid = None
        if resolver is not None:
            gid = resolver(amp)
        else:
            gid = (amp.get("glyma_id") or "").strip() or (amp.get("gene_id") or "").strip()

        if gid:
            if gid not in gene_ids:
                gene_ids.append(gid)
        else:
            missing.append(
                {
                    "subject_id": amp.get("subject_id", ""),
                    "gene_id": amp.get("gene_id", ""),
                    "glyma_id": amp.get("glyma_id", ""),
                }
            )

    return {
        "gene_ids": gene_ids,
        "missing_amplicon_gene_mapping": missing,
    }


def evaluate_blast_results(
    blast_results: List[Dict[str, Any]],
    expression_client: Any,
    resolver: Optional[Callable[[Dict[str, Any]], Optional[str]]] = None,
    use_cache: bool = True,
    attach_key: str = "e_eval",
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []

    for row in blast_results:
        item = dict(row)

        target_amplicons = row.get("target_amplicons", []) or []
        off_targets = row.get("off_targets", {}) or {}
        near_amplicons = off_targets.get("near", []) or []
        far_amplicons = off_targets.get("far", []) or []

        target_res = _resolve_amplicon_gene_ids(target_amplicons, resolver=resolver)
        near_res = _resolve_amplicon_gene_ids(near_amplicons, resolver=resolver)
        far_res = _resolve_amplicon_gene_ids(far_amplicons, resolver=resolver)

        target_expr = _merge_expr_for_gene_ids(
            target_res["gene_ids"],
            expression_client=expression_client,
            use_cache=use_cache,
        )
        near_expr = _merge_expr_for_gene_ids(
            near_res["gene_ids"],
            expression_client=expression_client,
            use_cache=use_cache,
        )
        far_expr = _merge_expr_for_gene_ids(
            far_res["gene_ids"],
            expression_client=expression_client,
            use_cache=use_cache,
        )

        target_expr_by_tissue = target_expr["expr_by_tissue"]
        near_expr_by_tissue = near_expr["expr_by_tissue"]
        far_expr_by_tissue = far_expr["expr_by_tissue"]

        target_expr_max = max(target_expr_by_tissue.values()) if target_expr_by_tissue else 0.0
        near_expr_max = max(near_expr_by_tissue.values()) if near_expr_by_tissue else 0.0
        far_expr_max = max(far_expr_by_tissue.values()) if far_expr_by_tissue else 0.0

        target_expr_sum = sum(target_expr_by_tissue.values()) if target_expr_by_tissue else 0.0
        near_expr_sum = sum(near_expr_by_tissue.values()) if near_expr_by_tissue else 0.0
        far_expr_sum = sum(far_expr_by_tissue.values()) if far_expr_by_tissue else 0.0

        # 展示规则优化：
        # 若 near/far 为 0，则 ratio 记为 None，而不是 8e8 那种不友好的数
        target_vs_near_ratio = None if near_expr_max <= 0 else target_expr_max / max(near_expr_max, EPS)
        target_vs_far_ratio = None if far_expr_max <= 0 else target_expr_max / max(far_expr_max, EPS)

        e_eval = {
            "target_genes": target_res["gene_ids"],
            "near_genes": near_res["gene_ids"],
            "far_genes": far_res["gene_ids"],

            "target_expr_by_tissue": target_expr_by_tissue,
            "near_expr_by_tissue": near_expr_by_tissue,
            "far_expr_by_tissue": far_expr_by_tissue,

            "target_expr_max": target_expr_max,
            "near_expr_max": near_expr_max,
            "far_expr_max": far_expr_max,

            "target_expr_sum": target_expr_sum,
            "near_expr_sum": near_expr_sum,
            "far_expr_sum": far_expr_sum,

            "dominant_target_tissues": _top_tissues(target_expr_by_tissue),
            "dominant_near_tissues": _top_tissues(near_expr_by_tissue),
            "dominant_far_tissues": _top_tissues(far_expr_by_tissue),

            "target_vs_near_ratio": target_vs_near_ratio,
            "target_vs_far_ratio": target_vs_far_ratio,

            "raw_detail": {
                "target": target_expr["raw_detail"],
                "near": near_expr["raw_detail"],
                "far": far_expr["raw_detail"],
            },
            "raw_to_canonical_detail": {
                "target": target_expr["raw_to_canonical_detail"],
                "near": near_expr["raw_to_canonical_detail"],
                "far": far_expr["raw_to_canonical_detail"],
            },

            "missing_expression_genes": {
                "target": target_expr["missing_expression_genes"],
                "near": near_expr["missing_expression_genes"],
                "far": far_expr["missing_expression_genes"],
            },
            "missing_amplicon_gene_mapping": {
                "target": target_res["missing_amplicon_gene_mapping"],
                "near": near_res["missing_amplicon_gene_mapping"],
                "far": far_res["missing_amplicon_gene_mapping"],
            },

            # 新增：交接/debug 时非常有用
            "expression_source_detail": {
                "target": target_expr["source_details"],
                "near": near_expr["source_details"],
                "far": far_expr["source_details"],
            },
        }

        item[attach_key] = e_eval
        out.append(item)

    return out