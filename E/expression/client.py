from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

from E.expression.cache import ExpressionCache
from E.expression.canonical import ConditionCanonicalizer
from E.expression.conekt_client import ConektExpressionClient, ensure_prefixed_glyma
from E.expression.soyomics_client import SoyOmicsExpressionClient
from E.expression.integrator import (
    SourceExpression,
    integrate_source_expressions,
)


@dataclass
class GeneExpression:
    glyma_gene: str
    prefixed_gene: str
    source: str
    ok: bool

    raw_expr: Dict[str, float] = field(default_factory=dict)
    raw_to_canonical: Dict[str, Optional[str]] = field(default_factory=dict)
    raw_canonical_details: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    canonical_expr: Dict[str, float] = field(default_factory=dict)

    notes: str = ""
    source_url: Optional[str] = None
    source_detail: Optional[str] = None

    # 新增：保留各 source 的输出，便于 debug
    by_source: Dict[str, Dict[str, Any]] = field(default_factory=dict)


class ExpressionClient:
    """
    Unified multi-source expression client.

    支持：
    - soyomics（当前拿到的是 macro-summary-like 数据）
    - conekt_lis
    - integrated（把多个 source 整合成一个统一输出）

    推荐：
    - preferred_source="integrated"
    """

    def __init__(
        self,
        cache_db_path: str = "expression_cache.sqlite",
        ttl_seconds: int = 7 * 24 * 3600,
        allowed_canon_conditions: Optional[Set[str]] = None,
        conekt_prefix: str = "glyma.Wm82.gnm2.ann1.",
        conekt_sleep_s: float = 0.2,
        conekt_timeout: int = 25,
        preferred_source: str = "integrated",
        enable_soyomics_fallback: bool = True,
        soyomics_timeout: int = 25,
        source_weights: Optional[Dict[str, float]] = None,
    ):
        self.cache = ExpressionCache(cache_db_path)
        self.ttl_seconds = ttl_seconds
        self.allowed = allowed_canon_conditions
        self.prefix = conekt_prefix
        self.preferred_source = preferred_source.strip().lower()
        self.enable_soyomics_fallback = enable_soyomics_fallback
        self.source_weights = source_weights or {}

        self.canon = ConditionCanonicalizer()
        self.conekt = ConektExpressionClient(
            timeout=conekt_timeout,
            sleep_s=conekt_sleep_s,
        )
        self.soyomics = SoyOmicsExpressionClient(timeout=soyomics_timeout)

    @staticmethod
    def _aggregate_canonical(values: List[float]) -> float:
        return max(values) if values else 0.0

    def _canonicalize_raw_expr(
        self,
        raw_expr: Dict[str, float],
        source: Optional[str],
    ) -> tuple[
        Dict[str, Optional[str]],
        Dict[str, Dict[str, Any]],
        Dict[str, float],
    ]:
        raw_to_canonical: Dict[str, Optional[str]] = {}
        raw_canonical_details: Dict[str, Dict[str, Any]] = {}
        canon_buckets: Dict[str, List[float]] = {}

        for raw_cond, value in raw_expr.items():
            detail = self.canon(raw_cond, source=source)
            canonical = detail["canonical_tissue"]

            if self.allowed is not None and canonical not in self.allowed:
                canonical_for_use: Optional[str] = None
            else:
                canonical_for_use = canonical

            raw_to_canonical[raw_cond] = canonical_for_use
            raw_canonical_details[raw_cond] = detail

            if canonical_for_use is not None:
                canon_buckets.setdefault(canonical_for_use, []).append(float(value))

        canonical_expr: Dict[str, float] = {
            tissue: self._aggregate_canonical(vals)
            for tissue, vals in canon_buckets.items()
        }

        return raw_to_canonical, raw_canonical_details, canonical_expr

    def _cache_source_key(self) -> str:
        return f"EXPR::{self.preferred_source}::soyfb{int(bool(self.enable_soyomics_fallback))}"

    def _from_cached(
        self,
        glyma_gene: str,
        prefixed: str,
        source: str,
        cached: Dict[str, Any],
    ) -> GeneExpression:
        raw_expr = {k: float(v) for k, v in (cached.get("raw_expr") or {}).items()}
        raw_to_canonical = dict(cached.get("raw_to_canonical") or {})
        raw_canonical_details = dict(cached.get("raw_canonical_details") or {})
        canonical_expr = {k: float(v) for k, v in (cached.get("canonical_expr") or {}).items()}

        return GeneExpression(
            glyma_gene=glyma_gene,
            prefixed_gene=prefixed,
            source=str(cached.get("source", source)),
            ok=bool(cached.get("ok", False)),
            raw_expr=raw_expr,
            raw_to_canonical=raw_to_canonical,
            raw_canonical_details=raw_canonical_details,
            canonical_expr=canonical_expr,
            notes=str(cached.get("notes", "cache_hit")),
            source_url=cached.get("source_url"),
            source_detail=cached.get("source_detail"),
            by_source=dict(cached.get("by_source") or {}),
        )

    def _build_gene_expression(
        self,
        glyma_gene: str,
        prefixed: str,
        source_name: str,
        res: Any,
    ) -> GeneExpression:
        raw_expr: Dict[str, float] = {}
        if res.ok and res.expr_by_condition:
            raw_expr = {k: float(v) for k, v in res.expr_by_condition.items()}

        raw_to_canonical, raw_canonical_details, canonical_expr = self._canonicalize_raw_expr(
            raw_expr,
            source=source_name,
        )

        return GeneExpression(
            glyma_gene=glyma_gene,
            prefixed_gene=prefixed,
            source=source_name,
            ok=res.ok,
            raw_expr=raw_expr,
            raw_to_canonical=raw_to_canonical,
            raw_canonical_details=raw_canonical_details,
            canonical_expr=canonical_expr,
            notes=res.notes,
            source_url=res.source_url,
            source_detail=source_name.upper(),
            by_source={
                source_name: {
                    "ok": res.ok,
                    "notes": res.notes,
                    "source_url": res.source_url,
                    "raw_expr": raw_expr,
                    "canonical_expr": canonical_expr,
                }
            },
        )

    def _fetch_from_conekt(self, glyma_gene: str) -> GeneExpression:
        prefixed = ensure_prefixed_glyma(glyma_gene, prefix=self.prefix)
        res = self.conekt.fetch_gene_expression(prefixed)
        return self._build_gene_expression(
            glyma_gene=glyma_gene,
            prefixed=prefixed,
            source_name="conekt_lis",
            res=res,
        )

    def _fetch_from_soyomics(self, glyma_gene: str) -> GeneExpression:
        prefixed = ensure_prefixed_glyma(glyma_gene, prefix=self.prefix)
        res = self.soyomics.fetch_gene_expression(glyma_gene)
        return self._build_gene_expression(
            glyma_gene=glyma_gene,
            prefixed=prefixed,
            source_name="soyomics",
            res=res,
        )

    def _integrate_gene(self, glyma_gene: str) -> GeneExpression:
        prefixed = ensure_prefixed_glyma(glyma_gene, prefix=self.prefix)

        soy = self._fetch_from_soyomics(glyma_gene)
        conekt = self._fetch_from_conekt(glyma_gene)

        integrated = integrate_source_expressions(
            [
                SourceExpression(
                    source="soyomics",
                    ok=soy.ok,
                    canonical_expr=soy.canonical_expr,
                    raw_expr=soy.raw_expr,
                    notes=soy.notes,
                    source_url=soy.source_url,
                ),
                SourceExpression(
                    source="conekt_lis",
                    ok=conekt.ok,
                    canonical_expr=conekt.canonical_expr,
                    raw_expr=conekt.raw_expr,
                    notes=conekt.notes,
                    source_url=conekt.source_url,
                ),
            ],
            source_weights=self.source_weights,
        )

        raw_expr: Dict[str, float] = {}
        raw_to_canonical: Dict[str, Optional[str]] = {}
        raw_canonical_details: Dict[str, Dict[str, Any]] = {}

        # integrated source 本身不保留混合后的 raw_expr，防止混淆
        by_source = {}
        by_source.update(soy.by_source)
        by_source.update(conekt.by_source)

        return GeneExpression(
            glyma_gene=glyma_gene,
            prefixed_gene=prefixed,
            source="integrated",
            ok=integrated.ok,
            raw_expr=raw_expr,
            raw_to_canonical=raw_to_canonical,
            raw_canonical_details=raw_canonical_details,
            canonical_expr=integrated.integrated_expr,
            notes=integrated.notes,
            source_url=None,
            source_detail="INTEGRATED",
            by_source=by_source,
        )

    def fetch_gene(self, glyma_gene: str, use_cache: bool = True) -> GeneExpression:
        """
        输入统一用 Glyma ID：
            Glyma.09G235300
        """
        prefixed = ensure_prefixed_glyma(glyma_gene, prefix=self.prefix)
        cache_source = self._cache_source_key()

        if use_cache:
            cached = self.cache.get(cache_source, prefixed, ttl_seconds=self.ttl_seconds)
            if cached:
                return self._from_cached(glyma_gene, prefixed, cache_source, cached)

        if self.preferred_source == "soyomics":
            ge = self._fetch_from_soyomics(glyma_gene)
            if not ge.ok and self.enable_soyomics_fallback:
                alt = self._fetch_from_conekt(glyma_gene)
                ge = alt if alt.ok else ge

        elif self.preferred_source == "conekt_lis":
            ge = self._fetch_from_conekt(glyma_gene)
            if not ge.ok and self.enable_soyomics_fallback:
                alt = self._fetch_from_soyomics(glyma_gene)
                ge = alt if alt.ok else ge

        else:
            ge = self._integrate_gene(glyma_gene)

        if use_cache:
            self.cache.set(
                cache_source,
                prefixed,
                {
                    "source": ge.source,
                    "ok": ge.ok,
                    "raw_expr": ge.raw_expr,
                    "raw_to_canonical": ge.raw_to_canonical,
                    "raw_canonical_details": ge.raw_canonical_details,
                    "canonical_expr": ge.canonical_expr,
                    "notes": ge.notes,
                    "source_url": ge.source_url,
                    "source_detail": ge.source_detail,
                    "by_source": ge.by_source,
                },
            )

        return ge

    def fetch_many(self, glyma_genes: List[str], use_cache: bool = True) -> List[GeneExpression]:
        seen = set()
        uniq: List[str] = []

        for g in glyma_genes:
            g0 = (g or "").strip()
            if g0 and g0 not in seen:
                seen.add(g0)
                uniq.append(g0)

        out: List[GeneExpression] = []
        for g in uniq:
            out.append(self.fetch_gene(g, use_cache=use_cache))
        return out

    def get_expr(self, glyma_gene: str, canonical_condition: str, use_cache: bool = True) -> float:
        ge = self.fetch_gene(glyma_gene, use_cache=use_cache)
        return float(ge.canonical_expr.get(canonical_condition, 0.0))