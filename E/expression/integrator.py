from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class SourceExpression:
    source: str
    ok: bool
    canonical_expr: Dict[str, float] = field(default_factory=dict)
    raw_expr: Dict[str, float] = field(default_factory=dict)
    notes: str = ""
    source_url: Optional[str] = None


@dataclass
class IntegratedExpression:
    ok: bool
    integrated_expr: Dict[str, float]
    by_source: Dict[str, SourceExpression]
    notes: str


DEFAULT_SOURCE_WEIGHTS: Dict[str, float] = {
    "soyomics": 1.0,
    "conekt_lis": 1.0,
}


def _max_merge(values: List[float]) -> float:
    return max(values) if values else 0.0


def integrate_source_expressions(
    expressions: List[SourceExpression],
    source_weights: Optional[Dict[str, float]] = None,
) -> IntegratedExpression:
    """
    把多个来源的 canonical expression 整合到同一个 tissue 空间。

    当前策略：
    1. 不同 source 各自先 canonicalize
    2. 对同一个 tissue 采用 weighted max
    3. 保留每个 source 的原始结果，供 debug / 交接

    说明：
    - 这里先用 weighted max，而不是平均值
    - 原因是不同数据库量纲与处理流程未必完全一致
    - weighted max 更适合“是否支持该组织表达”这个场景
    """
    weights = dict(DEFAULT_SOURCE_WEIGHTS)
    if source_weights:
        weights.update(source_weights)

    by_source: Dict[str, SourceExpression] = {}
    tissue_buckets: Dict[str, List[float]] = {}

    ok_any = False
    notes_list: List[str] = []

    for expr in expressions:
        by_source[expr.source] = expr
        notes_list.append(f"{expr.source}: {expr.notes}")

        if not expr.ok:
            continue

        ok_any = True
        w = float(weights.get(expr.source, 1.0))

        for tissue, value in expr.canonical_expr.items():
            tissue_buckets.setdefault(tissue, []).append(float(value) * w)

    integrated_expr: Dict[str, float] = {
        tissue: _max_merge(vals)
        for tissue, vals in tissue_buckets.items()
    }

    return IntegratedExpression(
        ok=ok_any,
        integrated_expr=integrated_expr,
        by_source=by_source,
        notes=" | ".join(notes_list),
    )