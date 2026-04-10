from __future__ import annotations

"""
canonical.py
============

本模块负责把表达数据库中的 raw condition 标签，解析成项目内部统一使用的
canonical tissue，并保留原始标签与元信息，供展示、追溯和 debug 使用。

【设计原则 / 交接警告】
1. 采用“双层表示”原则：
   - raw condition 必须完整保留，用于追溯、展示和调试；
   - canonical tissue 用于表达量聚合、组织比较与最终打分。
   最终排名只基于 canonical tissue，不直接基于 raw condition。

2. 不丢信息：
   即使 condition 无法识别，也必须保留：
   - raw
   - tokens
   - unmatched_tokens
   并把 canonical_tissue 设为 "unknown"。

3. 不要强行猜：
   如果不能可靠判断，就返回 unknown，而不是硬猜成某个组织。

4. 当前支持两套主要表达来源：
   - CoNekT/LIS 老 atlas（root hair / Stacey）
   - SoyOmics 28 tissues

5. source-aware 规则：
   某些标签（比如 stem / leafbud）只有在 SoyOmics 28 tissues 里才有明确组织学含义。
"""

import re
from typing import Any, Dict, List, Optional, Set


CANONICAL_TISSUES: Set[str] = {
    # 老 atlas
    "root_hair",
    "root_tip",
    "root",
    "nodule",
    "leaf",
    "flower",
    "pod",
    "sam",
    # SoyOmics 28 tissues 扩展
    "cotyledon",
    "hypocotyl",
    "shoot_apex",
    "simple_leaf",
    "trifoliate_leaf",
    "pod_seed",
    "seed",
    "axillary_bud",
    # fallback
    "unknown",
}

TISSUE_PRIORITY: List[str] = [
    "root_hair",
    "root_tip",
    "nodule",
    "sam",
    "shoot_apex",
    "axillary_bud",
    "cotyledon",
    "hypocotyl",
    "simple_leaf",
    "trifoliate_leaf",
    "pod_seed",
    "seed",
    "root",
    "leaf",
    "flower",
    "pod",
]

TOKEN_TO_TISSUE_BASE: Dict[str, str] = {
    # root hair
    "RH": "root_hair",
    "ROOTHAIR": "root_hair",
    "ROOT_HAIR": "root_hair",

    # root tip
    "ROOTTIP": "root_tip",
    "ROOT_TIP": "root_tip",

    # root
    "ROOT": "root",
    "ROOTS": "root",
    "STRIPPEDROOT": "root",
    "STRIPPED_ROOT": "root",

    # nodule
    "NODULE": "nodule",
    "NODULES": "nodule",
    "NOD": "nodule",

    # leaf
    "LEAF": "leaf",
    "LEAVES": "leaf",

    # flower
    "FLOWER": "flower",
    "FLOWERS": "flower",
    "FLORAL": "flower",

    # pod
    "POD": "pod",
    "PODS": "pod",
    "GREENPOD": "pod",
    "GREEN_POD": "pod",

    # shoot apical meristem / sam
    "SAM": "sam",
    "MERISTEM": "sam",
    "SHOOTAPICALMERISTEM": "sam",
    "SHOOT_APICAL_MERISTEM": "sam",
    "APICALMERISTEM": "sam",
}

TOKEN_TO_TISSUE_SOYOMICS: Dict[str, str] = {
    # SoyOmics 28 tissues
    "COTYLEDON": "cotyledon",
    "SEEDLINGCOTYLEDON": "cotyledon",
    "HYPOCOTYL": "hypocotyl",
    "STEM": "hypocotyl",          # source-specific: SoyOmics 中 stem 归到下胚轴/茎段体系
    "LEAFBUD": "shoot_apex",      # source-specific
    "SHOOTAPEX": "shoot_apex",
    "SHOOT_APEX": "shoot_apex",
    "AXILLARYBUD": "axillary_bud",
    "AXILLARY_BUD": "axillary_bud",
    "SEED": "seed",
    "SEEDS": "seed",
    "PODSEED": "pod_seed",
    "POD_SEED": "pod_seed",
    "POD&SEED": "pod_seed",
    "SIMPLELEAF": "simple_leaf",
    "SIMPLE_LEAF": "simple_leaf",
    "TRIFOLIATELEAF": "trifoliate_leaf",
    "TRIFOLIATE_LEAF": "trifoliate_leaf",
}

INFECTION_TOKENS: Set[str] = {
    "IN",
    "INFECTED",
    "INOCULATED",
    "INFECTION",
}

MOCK_TOKENS: Set[str] = {
    "MOCK",
    "CONTROL",
    "CTRL",
    "UN",
    "UNINFECTED",
    "UNTREATED",
}

NEUTRAL_TOKENS: Set[str] = {
    "CELL",
    "CELLS",
    "MATURE",
    "GREEN",
    "SHOOT",
    "APICAL",
    "HAIR",
    "TIP",
    "STRIPPED",
    "STACEY",
    "SCRIP",
    "BUD",        # 单独出现时不作为噪声报错；真正组织判断靠组合规则
    "AXILLARY",   # 同上
}

TIMEPOINT_RE = re.compile(r"^(?P<time>\d+)(?P<unit>HAS|HAI)$", re.IGNORECASE)
COMPACT_HA_RE = re.compile(r"^(?P<time>\d+)HA(?P<rep>\d+)$", re.IGNORECASE)
REPLICATE_PATTERNS = [
    re.compile(r"^HA(?P<rep>\d+)$", re.IGNORECASE),
    re.compile(r"^REP(?P<rep>\d+)$", re.IGNORECASE),
    re.compile(r"^R(?P<rep>\d+)$", re.IGNORECASE),
    re.compile(r"^(?P<rep>\d+)$", re.IGNORECASE),  # soyomics labels like stem-1
]


def _normalize_raw(raw: str) -> str:
    return raw.strip()


def _compact_token(token: str) -> str:
    return re.sub(r"[^A-Z0-9&]", "", token.upper())


def tokenize_condition(raw: str) -> List[str]:
    if not isinstance(raw, str):
        raw = str(raw)
    raw = _normalize_raw(raw).upper()
    parts = re.split(r"[\s_\-\/\(\)\[\],;:]+", raw)
    return [p for p in parts if p]


def _source_key(source: Optional[str]) -> Optional[str]:
    if source is None:
        return None
    s = str(source).strip().lower()
    if s in {"soyomics", "soyomics_28tissues", "soyomics_28_tissues"}:
        return "soyomics"
    if s in {"conekt", "conekt_lis", "lis", "conekt/lis"}:
        return "conekt_lis"
    return s


def _get_token_map(source: Optional[str]) -> Dict[str, str]:
    sk = _source_key(source)
    if sk == "soyomics":
        merged = dict(TOKEN_TO_TISSUE_BASE)
        merged.update(TOKEN_TO_TISSUE_SOYOMICS)
        return merged
    return dict(TOKEN_TO_TISSUE_BASE)


def extract_metadata(tokens: List[str]) -> Dict[str, Any]:
    metadata: Dict[str, Any] = {
        "timepoint": None,
        "time_unit": None,
        "infection": None,
        "mock": None,
        "replicate": None,
        "notes": [],
    }

    for tok in tokens:
        tok_u = tok.upper()
        tok_c = _compact_token(tok_u)

        m = TIMEPOINT_RE.match(tok_c)
        if m:
            metadata["timepoint"] = m.group("time")
            metadata["time_unit"] = m.group("unit").upper()
            continue

        m = COMPACT_HA_RE.match(tok_c)
        if m:
            if metadata["timepoint"] is None:
                metadata["timepoint"] = m.group("time")
            if metadata["time_unit"] is None:
                metadata["time_unit"] = "HAI"
            if metadata["replicate"] is None:
                metadata["replicate"] = m.group("rep")
            continue

        for pat in REPLICATE_PATTERNS:
            mm = pat.match(tok_c)
            if mm:
                # 避免把完整 timepoint 误认成 replicate
                if tok_c.isdigit() and len(tok_c) >= 2 and metadata["timepoint"] is None:
                    break
                if metadata["replicate"] is None:
                    metadata["replicate"] = mm.group("rep")
                break

        if tok_c in INFECTION_TOKENS:
            metadata["infection"] = True
            if metadata["mock"] is None:
                metadata["mock"] = False
            continue

        if tok_c in MOCK_TOKENS:
            metadata["mock"] = True
            if metadata["infection"] is None:
                metadata["infection"] = False
            continue

    return metadata


def _add_hit(hits: Dict[str, List[str]], tissue: str, token: str) -> None:
    hits.setdefault(tissue, [])
    if token not in hits[tissue]:
        hits[tissue].append(token)


def _match_tissues(tokens: List[str], source: Optional[str] = None) -> Dict[str, List[str]]:
    hits: Dict[str, List[str]] = {}
    token_map = _get_token_map(source)

    token_set = {t.upper() for t in tokens}
    compact_set = {_compact_token(t) for t in tokens}
    sk = _source_key(source)

    # -------------------------
    # base combination rules
    # -------------------------
    if ("ROOT" in token_set and "HAIR" in token_set) or ("ROOT" in compact_set and "HAIR" in compact_set):
        _add_hit(hits, "root_hair", "ROOT")
        _add_hit(hits, "root_hair", "HAIR")

    if ("ROOT" in token_set and "TIP" in token_set) or ("ROOT" in compact_set and "TIP" in compact_set):
        _add_hit(hits, "root_tip", "ROOT")
        _add_hit(hits, "root_tip", "TIP")

    has_stripped = "STRIPPED" in token_set or "STRIPPED" in compact_set
    has_root = "ROOT" in token_set or "ROOTS" in token_set or "ROOT" in compact_set or "ROOTS" in compact_set
    if has_stripped and has_root:
        _add_hit(hits, "root", "STRIPPED")
        _add_hit(hits, "root", "ROOT")

    has_green = "GREEN" in token_set or "GREEN" in compact_set
    has_pod = "POD" in token_set or "PODS" in token_set or "POD" in compact_set or "PODS" in compact_set
    if has_green and has_pod:
        _add_hit(hits, "pod", "GREEN")
        _add_hit(hits, "pod", "POD")

    if (
        ("SHOOT" in token_set or "SHOOT" in compact_set)
        and ("APICAL" in token_set or "APICAL" in compact_set)
        and ("MERISTEM" in token_set or "MERISTEM" in compact_set)
    ):
        _add_hit(hits, "sam", "SHOOT")
        _add_hit(hits, "sam", "APICAL")
        _add_hit(hits, "sam", "MERISTEM")

    # -------------------------
    # SoyOmics-specific combination rules
    # -------------------------
    if sk == "soyomics":
        # pod + seed -> pod_seed
        if ("POD" in token_set and "SEED" in token_set) or ("POD&SEED" in compact_set):
            _add_hit(hits, "pod_seed", "POD")
            _add_hit(hits, "pod_seed", "SEED")

        # axillary + bud -> axillary_bud
        if ("AXILLARY" in token_set and "BUD" in token_set) or ("AXILLARYBUD" in compact_set):
            _add_hit(hits, "axillary_bud", "AXILLARY")
            _add_hit(hits, "axillary_bud", "BUD")

    # -------------------------
    # single token hits
    # -------------------------
    for tok in tokens:
        tok_u = tok.upper()
        tok_c = _compact_token(tok_u)

        if tok_u in token_map:
            _add_hit(hits, token_map[tok_u], tok_u)
        if tok_c in token_map:
            _add_hit(hits, token_map[tok_c], tok_u)

    return hits


def _choose_canonical_tissue(tissue_hits: Dict[str, List[str]]) -> str:
    if not tissue_hits:
        return "unknown"
    for tissue in TISSUE_PRIORITY:
        if tissue in tissue_hits:
            return tissue
    return "unknown"


def _estimate_confidence(canonical_tissue: str, matched_tokens: List[str], unmatched_tokens: List[str]) -> str:
    if canonical_tissue == "unknown":
        return "low"
    if matched_tokens and len(unmatched_tokens) <= 2:
        return "high"
    return "medium"


def canonicalize_condition(raw: str, source: Optional[str] = None) -> Dict[str, Any]:
    raw_norm = _normalize_raw(raw)
    tokens = tokenize_condition(raw_norm)
    metadata = extract_metadata(tokens)

    tissue_hits = _match_tissues(tokens, source=source)
    canonical_tissue = _choose_canonical_tissue(tissue_hits)

    matched_tokens: List[str] = []
    if canonical_tissue in tissue_hits:
        matched_tokens = tissue_hits[canonical_tissue][:]

    token_map = _get_token_map(source)
    matched_upper = {t.upper() for t in matched_tokens}
    unmatched_tokens: List[str] = []

    for tok in tokens:
        tok_u = tok.upper()
        tok_c = _compact_token(tok_u)

        if tok_u in matched_upper:
            continue

        if tok_u in token_map or tok_c in token_map:
            continue

        if TIMEPOINT_RE.match(tok_c):
            continue
        if COMPACT_HA_RE.match(tok_c):
            continue
        if tok_c in INFECTION_TOKENS or tok_c in MOCK_TOKENS or tok_c in NEUTRAL_TOKENS:
            continue
        if any(p.match(tok_c) for p in REPLICATE_PATTERNS):
            continue

        unmatched_tokens.append(tok)

    if canonical_tissue == "unknown":
        metadata["notes"].append("unrecognized condition")

    confidence = _estimate_confidence(
        canonical_tissue=canonical_tissue,
        matched_tokens=matched_tokens,
        unmatched_tokens=unmatched_tokens,
    )

    return {
        "raw": raw_norm,
        "tokens": tokens,
        "canonical_tissue": canonical_tissue,
        "matched_tokens": matched_tokens,
        "unmatched_tokens": unmatched_tokens,
        "metadata": metadata,
        "confidence": confidence,
        "source": _source_key(source),
    }


def canonical_condition(raw: str, source: Optional[str] = None) -> str:
    return canonicalize_condition(raw, source=source)["canonical_tissue"]


class ConditionCanonicalizer:
    def __call__(self, raw: str, source: Optional[str] = None) -> Dict[str, Any]:
        return canonicalize_condition(raw, source=source)

    def canonical(self, raw: str, source: Optional[str] = None) -> str:
        return canonical_condition(raw, source=source)

    def batch(self, raws: List[str], source: Optional[str] = None) -> List[Dict[str, Any]]:
        return [canonicalize_condition(r, source=source) for r in raws]


if __name__ == "__main__":
    import json

    examples = [
        # old atlas
        "12HA1_IN_RH",
        "84HAS RH",
        "120HAS root hair",
        "root tip",
        "ROOT",
        "mature nodules",
        "leaf",
        "flower",
        "green pod",
        "shoot apical meristem",
        "24HAI mock RH",
        "48HAI inoculated stripped roots",
        # soyomics
        "cotyledon-1",
        "stem-1",
        "leafbud-2",
        "leaf-3",
        "flower-5",
        "pod-3",
        "seed-2",
        "axillary_bud",
        "axillary bud",
        "root",
        "XYZ_ABC",
    ]

    for ex in examples:
        print("=" * 80)
        print(json.dumps(canonicalize_condition(ex, source="soyomics"), ensure_ascii=False, indent=2))