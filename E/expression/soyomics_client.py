from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import requests


BASE_URL = "https://ngdc.cncb.ac.cn/soyomics/transcriptome"

HEADERS = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "Referer": "https://ngdc.cncb.ac.cn/soyomics/transcriptome/tissues",
    "X-Requested-With": "XMLHttpRequest",
}

# 页面上可见的 28 samples 标签
KNOWN_28_SAMPLES: List[str] = [
    "cotyledon-1",
    "cotyledon-2",
    "stem-1",
    "stem-2",
    "leafbud-1",
    "leafbud-2",
    "leafbud-3",
    "leaf-1",
    "leaf-2",
    "leaf-3",
    "flo-1",
    "flo-2",
    "flo-3",
    "flo-4",
    "flo-5",
    "pod & seed-1",
    "pod & seed-2",
    "pod & seed-3",
    "pod-1",
    "pod-2",
    "pod-3",
    "seed-1",
    "seed-2",
    "seed-3",
    "seed-4",
    "seed-5",
    "axillary_bud",
    "root",
]

# 从 other_List URL / 返回结构中能推断出来的 summary 级列
KNOWN_MACRO_TISSUES: List[str] = [
    "stem",
    "young_leaf",
    "mature_leaf",
    "old_leaf",
    "flower",
    "pod_seed",
    "seed6",
    "seed8",
    "axillary_bud",
    "root",
]

KNOWN_28_SAMPLES_SET = {x.lower() for x in KNOWN_28_SAMPLES}
KNOWN_MACRO_TISSUES_SET = {x.lower() for x in KNOWN_MACRO_TISSUES}

IGNORE_METADATA_KEYS = {
    "gene_id",
    "id",
    "name",
    "zh13_id",
    "wm82_id",
    "locus",
    "gene_start",
    "gene_end",
    "chr",
    "chromosome",
    "strand",
}


@dataclass
class SoyOmicsExpressionResult:
    ok: bool
    expr_by_condition: Dict[str, float]
    notes: str
    source_url: Optional[str] = None


def _is_number(x: Any) -> bool:
    try:
        float(x)
        return True
    except Exception:
        return False


def _normalize_label(x: Any) -> str:
    return str(x).strip()


def _safe_preview(x: Any, max_items: int = 3) -> str:
    try:
        if isinstance(x, list):
            return repr(x[:max_items])
        if isinstance(x, dict):
            ks = list(x.keys())[:max_items]
            return repr({k: x[k] for k in ks})
        return repr(x)
    except Exception:
        return "<preview_error>"


def _recognized_key_score(keys: List[str]) -> tuple[int, int]:
    """
    返回：
    - 28 sample 命中数
    - macro tissue 命中数
    """
    s28 = sum(1 for k in keys if k.lower() in KNOWN_28_SAMPLES_SET)
    smacro = sum(1 for k in keys if k.lower() in KNOWN_MACRO_TISSUES_SET)
    return s28, smacro


def _parse_categories_series_dict(data: Dict[str, Any]) -> Dict[str, float]:
    """
    解析：
    {
        "categories": [...],
        "series": [{"data": [...]}]
    }
    """
    expr: Dict[str, float] = {}

    categories = data.get("categories")
    series = data.get("series")
    if not (isinstance(categories, list) and isinstance(series, list) and series):
        return expr

    first = series[0]
    if isinstance(first, dict):
        values = first.get("data")
        if isinstance(values, list):
            for tissue, value in zip(categories, values):
                if _is_number(value):
                    expr[_normalize_label(tissue)] = float(value)

    return expr


def _parse_samples_data_dict(data: Dict[str, Any]) -> Dict[str, float]:
    """
    解析：
    {
        "samples": [...],
        "data": [...]
    }
    或
    {
        "labels": [...],
        "values": [...]
    }
    """
    expr: Dict[str, float] = {}

    pairs = [
        ("samples", "data"),
        ("labels", "values"),
        ("x", "y"),
    ]

    for lk, vk in pairs:
        labels = data.get(lk)
        values = data.get(vk)
        if isinstance(labels, list) and isinstance(values, list):
            tmp: Dict[str, float] = {}
            for tissue, value in zip(labels, values):
                if _is_number(value):
                    tmp[_normalize_label(tissue)] = float(value)
            if tmp:
                return tmp

    return expr


def _parse_row_dict_as_mapping(row: Dict[str, Any]) -> Dict[str, float]:
    """
    把一个“表格行 dict”解释成表达映射，但只接受：
    - 28 sample keys
    - macro tissue keys

    明确排除 gene_start / gene_end 等元数据字段。
    """
    expr: Dict[str, float] = {}

    for k, v in row.items():
        kl = k.lower().strip()
        if kl in IGNORE_METADATA_KEYS:
            continue
        if kl in KNOWN_28_SAMPLES_SET or kl in KNOWN_MACRO_TISSUES_SET:
            if _is_number(v):
                expr[_normalize_label(k)] = float(v)

    return expr


def _parse_list_of_pairs(data: List[Any]) -> Dict[str, float]:
    """
    解析：
    [
        ["stem", 12.3],
        ["young_leaf", 0.1],
        ...
    ]
    """
    expr: Dict[str, float] = {}

    if not data:
        return expr

    if all(isinstance(item, list) and len(item) >= 2 for item in data):
        for item in data:
            label = _normalize_label(item[0])
            value = item[1]
            if _is_number(value):
                expr[label] = float(value)

    return expr


def _parse_list_of_dicts_shallow(data: List[Any]) -> Dict[str, float]:
    """
    解析：
    [
        {"name":"stem","value":12.3},
        {"sample":"cotyledon-1","expr":0.4},
        ...
    ]
    """
    expr: Dict[str, float] = {}

    if not data or not all(isinstance(item, dict) for item in data):
        return expr

    label_keys = ("name", "sample", "label", "condition", "tissue", "x")
    value_keys = ("value", "expr", "expression", "y", "mean")

    for item in data:
        label = None
        value = None

        for lk in label_keys:
            if lk in item and str(item[lk]).strip():
                label = item[lk]
                break

        for vk in value_keys:
            if vk in item and _is_number(item[vk]):
                value = item[vk]
                break

        if label is not None and value is not None:
            expr[_normalize_label(label)] = float(value)

    return expr


def _parse_two_parallel_lists(data: List[Any]) -> Dict[str, float]:
    """
    解析：
    [
        [labels...],
        [values...]
    ]
    """
    expr: Dict[str, float] = {}

    if len(data) != 2:
        return expr
    if not isinstance(data[0], list) or not isinstance(data[1], list):
        return expr

    labels = data[0]
    values = data[1]

    for tissue, value in zip(labels, values):
        if _is_number(value):
            expr[_normalize_label(tissue)] = float(value)

    return expr


def _parse_numeric_vector_with_known_order(data: List[Any]) -> Dict[str, float]:
    """
    如果 top-level 是纯数值 list，则尝试用已知顺序解释。
    """
    expr: Dict[str, float] = {}

    if not data:
        return expr
    if not all(_is_number(x) for x in data):
        return expr

    if len(data) == len(KNOWN_28_SAMPLES):
        for tissue, value in zip(KNOWN_28_SAMPLES, data):
            expr[tissue] = float(value)
        return expr

    if len(data) == len(KNOWN_MACRO_TISSUES):
        for tissue, value in zip(KNOWN_MACRO_TISSUES, data):
            expr[tissue] = float(value)
        return expr

    return expr


def _prefer_expr(a: Dict[str, float], b: Dict[str, float]) -> Dict[str, float]:
    """
    选择更像“真正表达矩阵”的那个：
    优先 28-sample 命中数高的；
    若相同，选择条目数更多的。
    """
    if not a:
        return b
    if not b:
        return a

    a28, amacro = _recognized_key_score(list(a.keys()))
    b28, bmacro = _recognized_key_score(list(b.keys()))

    if b28 > a28:
        return b
    if a28 > b28:
        return a

    if len(b) > len(a):
        return b
    return a


def _parse_any(data: Any) -> Dict[str, float]:
    """
    递归解析，优先寻找更深层的 28-sample 数据，
    最后才退回到 summary row。
    """
    best: Dict[str, float] = {}

    if isinstance(data, dict):
        # 1. 先试显式 chart 结构
        expr = _parse_categories_series_dict(data)
        best = _prefer_expr(best, expr)

        expr = _parse_samples_data_dict(data)
        best = _prefer_expr(best, expr)

        # 2. 再递归深入 dict values，优先找更深层 28-sample 数据
        for v in data.values():
            expr = _parse_any(v)
            best = _prefer_expr(best, expr)

        # 3. 最后才把当前 dict 当 row dict
        expr = _parse_row_dict_as_mapping(data)
        best = _prefer_expr(best, expr)

        return best

    if isinstance(data, list):
        # 1. 先试标准 list 结构
        expr = _parse_two_parallel_lists(data)
        best = _prefer_expr(best, expr)

        expr = _parse_list_of_pairs(data)
        best = _prefer_expr(best, expr)

        expr = _parse_list_of_dicts_shallow(data)
        best = _prefer_expr(best, expr)

        expr = _parse_numeric_vector_with_known_order(data)
        best = _prefer_expr(best, expr)

        # 2. 再递归深入 list items
        for item in data:
            expr = _parse_any(item)
            best = _prefer_expr(best, expr)

        return best

    return {}


class SoyOmicsExpressionClient:
    """
    SoyOmics expression client

    已确认页面实际调用：
        /transcriptome/other_List/heatmap?gene_id=<GlymaID>
    """

    def __init__(self, timeout: int = 25):
        self.timeout = timeout

    def fetch_gene_expression(self, glyma_id: str) -> SoyOmicsExpressionResult:
        url = f"{BASE_URL}/other_List/heatmap"
        params = {"gene_id": glyma_id}

        try:
            r = requests.get(
                url,
                params=params,
                headers=HEADERS,
                timeout=self.timeout,
            )
        except Exception as e:
            return SoyOmicsExpressionResult(
                ok=False,
                expr_by_condition={},
                notes=f"request error: {e}",
                source_url=url,
            )

        if r.status_code != 200:
            return SoyOmicsExpressionResult(
                ok=False,
                expr_by_condition={},
                notes=f"HTTP {r.status_code}",
                source_url=r.url,
            )

        try:
            data = r.json()
        except Exception as e:
            return SoyOmicsExpressionResult(
                ok=False,
                expr_by_condition={},
                notes=f"json parse error: {e}",
                source_url=r.url,
            )

        expr = _parse_any(data)

        if not expr:
            top_type = type(data).__name__
            preview = _safe_preview(data)
            return SoyOmicsExpressionResult(
                ok=False,
                expr_by_condition={},
                notes=f"empty expression; top-level type={top_type}; preview={preview}",
                source_url=r.url,
            )

        keys = list(expr.keys())
        hit28, hitmacro = _recognized_key_score(keys)

        # 给 notes 明确标记：拿到的是 28-sample 还是 macro-summary
        if hit28 >= 6:
            note = f"ok (28-sample-like; matched_28={hit28})"
        elif hitmacro >= 3:
            note = f"ok (macro-summary-like; matched_macro={hitmacro})"
        else:
            note = f"ok (unclassified; matched_28={hit28}; matched_macro={hitmacro})"

        return SoyOmicsExpressionResult(
            ok=True,
            expr_by_condition=expr,
            notes=note,
            source_url=r.url,
        )