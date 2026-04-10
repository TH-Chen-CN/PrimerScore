from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]  # .../new
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from E.evaluate import (  # noqa: E402
    build_expression_summary_for_pair,
    evaluate_blast_results,
)


def _pretty(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True, default=str)


class FakeExpressionClient:
    """
    假 expression client：
    fetch_many(gene_ids) -> list[GeneExpression-like]
    """
    def __init__(self, gene_expr_map):
        self.gene_expr_map = gene_expr_map
        self.calls = []

    def fetch_many(self, gene_ids, use_cache=True):
        self.calls.append((list(gene_ids), use_cache))
        out = []
        for gid in gene_ids:
            payload = self.gene_expr_map.get(gid)
            if payload is None:
                out.append(
                    SimpleNamespace(
                        glyma_gene=gid,
                        ok=False,
                        raw_expr={},
                        raw_to_canonical={},
                        raw_canonical_details={},
                        canonical_expr={},
                        notes="missing",
                        source_url=None,
                    )
                )
            else:
                out.append(
                    SimpleNamespace(
                        glyma_gene=gid,
                        ok=True,
                        raw_expr=payload["raw_expr"],
                        raw_to_canonical=payload["raw_to_canonical"],
                        raw_canonical_details=payload.get("raw_canonical_details", {}),
                        canonical_expr=payload["canonical_expr"],
                        notes="ok",
                        source_url=None,
                    )
                )
        return out


def fake_resolver(amplicon):
    """
    优先用 glyma_id；否则按 subject_id 映射。
    """
    if amplicon.get("glyma_id"):
        return amplicon["glyma_id"]

    mapping = {
        "XM_TARGET": "Glyma.TARGET",
        "XM_NEAR": "Glyma.NEAR",
        "XM_FAR": "Glyma.FAR",
    }
    return mapping.get(amplicon.get("subject_id"))


def make_fake_pair():
    return {
        "pair_index": 0,
        "forward": "AAAA",
        "reverse": "TTTT",
        "product_size": 708,
        "tm_diff": 0.5,
        "heterodimer_dG": -2.1,
        "target_subject_id": "XM_TARGET",
        "target_amplicons": [
            {
                "subject_id": "XM_TARGET",
                "gene_id": "XM_TARGET",
                "glyma_id": "",
                "length": 708,
            }
        ],
        "off_targets": {
            "near": [
                {
                    "subject_id": "XM_NEAR",
                    "gene_id": "XM_NEAR",
                    "glyma_id": "",
                    "length": 720,
                }
            ],
            "far": [
                {
                    "subject_id": "XM_FAR",
                    "gene_id": "XM_FAR",
                    "glyma_id": "",
                    "length": 1100,
                }
            ],
        },
    }


def test_build_expression_summary_for_pair_basic():
    pair = make_fake_pair()

    client = FakeExpressionClient(
        {
            "Glyma.TARGET": {
                "raw_expr": {
                    "84HAS RH": 5.0,
                    "120HAS root hair": 7.0,
                    "ROOT": 3.2,
                },
                "raw_to_canonical": {
                    "84HAS RH": "root_hair",
                    "120HAS root hair": "root_hair",
                    "ROOT": "root",
                },
                "canonical_expr": {
                    "root_hair": 7.0,
                    "root": 3.2,
                },
            },
            "Glyma.NEAR": {
                "raw_expr": {
                    "ROOT": 2.1,
                    "LEAF": 0.4,
                },
                "raw_to_canonical": {
                    "ROOT": "root",
                    "LEAF": "leaf",
                },
                "canonical_expr": {
                    "root": 2.1,
                    "leaf": 0.4,
                },
            },
            "Glyma.FAR": {
                "raw_expr": {
                    "FLOWER": 0.8,
                },
                "raw_to_canonical": {
                    "FLOWER": "flower",
                },
                "canonical_expr": {
                    "flower": 0.8,
                },
            },
        }
    )

    e_eval = build_expression_summary_for_pair(
        pair,
        expression_client=client,
        resolver=fake_resolver,
        use_cache=False,
    )

    assert e_eval["target_genes"] == ["Glyma.TARGET"], _pretty(e_eval)
    assert e_eval["near_genes"] == ["Glyma.NEAR"], _pretty(e_eval)
    assert e_eval["far_genes"] == ["Glyma.FAR"], _pretty(e_eval)

    assert e_eval["target_expr_by_tissue"] == {"root_hair": 7.0, "root": 3.2}, _pretty(e_eval)
    assert e_eval["near_expr_by_tissue"] == {"root": 2.1, "leaf": 0.4}, _pretty(e_eval)
    assert e_eval["far_expr_by_tissue"] == {"flower": 0.8}, _pretty(e_eval)

    assert e_eval["target_expr_max"] == 7.0, _pretty(e_eval)
    assert e_eval["near_expr_max"] == 2.1, _pretty(e_eval)
    assert e_eval["far_expr_max"] == 0.8, _pretty(e_eval)

    assert e_eval["target_expr_sum"] == 10.2, _pretty(e_eval)
    assert e_eval["near_expr_sum"] == 2.5, _pretty(e_eval)
    assert e_eval["far_expr_sum"] == 0.8, _pretty(e_eval)

    assert abs(e_eval["target_vs_near_ratio"] - (7.0 / 2.1)) < 1e-9, _pretty(e_eval)
    assert abs(e_eval["target_vs_far_ratio"] - (7.0 / 0.8)) < 1e-9, _pretty(e_eval)

    assert e_eval["dominant_target_tissues"] == ["root_hair", "root"], _pretty(e_eval)
    assert e_eval["dominant_near_tissues"] == ["root", "leaf"], _pretty(e_eval)
    assert e_eval["dominant_far_tissues"] == ["flower"], _pretty(e_eval)

    assert "Glyma.TARGET" in e_eval["raw_detail"]["target"], _pretty(e_eval)
    assert "Glyma.NEAR" in e_eval["raw_detail"]["near"], _pretty(e_eval)
    assert "Glyma.FAR" in e_eval["raw_detail"]["far"], _pretty(e_eval)


def test_group_internal_same_tissue_uses_max():
    """
    同一组内，如果两个 transcript / gene 都落到同一个 tissue，
    evaluate 规则应使用 max，而不是 sum。
    """
    pair = {
        "target_amplicons": [
            {"subject_id": "XM_A", "gene_id": "XM_A", "glyma_id": "Glyma.A", "length": 700},
            {"subject_id": "XM_B", "gene_id": "XM_B", "glyma_id": "Glyma.B", "length": 700},
        ],
        "off_targets": {"near": [], "far": []},
    }

    client = FakeExpressionClient(
        {
            "Glyma.A": {
                "raw_expr": {"84HAS RH": 5.0},
                "raw_to_canonical": {"84HAS RH": "root_hair"},
                "canonical_expr": {"root_hair": 5.0},
            },
            "Glyma.B": {
                "raw_expr": {"120HAS root hair": 7.0},
                "raw_to_canonical": {"120HAS root hair": "root_hair"},
                "canonical_expr": {"root_hair": 7.0},
            },
        }
    )

    e_eval = build_expression_summary_for_pair(
        pair,
        expression_client=client,
        resolver=fake_resolver,
        use_cache=False,
    )

    assert e_eval["target_expr_by_tissue"] == {"root_hair": 7.0}, _pretty(e_eval)
    assert e_eval["target_expr_sum"] == 7.0, _pretty(e_eval)
    assert e_eval["target_expr_max"] == 7.0, _pretty(e_eval)


def test_missing_mapping_and_missing_expression_are_preserved():
    pair = {
        "target_amplicons": [
            {"subject_id": "XM_UNKNOWN", "gene_id": "XM_UNKNOWN", "glyma_id": "", "length": 700},
            {"subject_id": "XM_TARGET", "gene_id": "XM_TARGET", "glyma_id": "", "length": 700},
        ],
        "off_targets": {"near": [], "far": []},
    }

    client = FakeExpressionClient(
        {
            "Glyma.TARGET": {
                "raw_expr": {"ROOT": 3.0},
                "raw_to_canonical": {"ROOT": "root"},
                "canonical_expr": {"root": 3.0},
            }
        }
    )

    e_eval = build_expression_summary_for_pair(
        pair,
        expression_client=client,
        resolver=fake_resolver,
        use_cache=False,
    )

    assert e_eval["target_genes"] == ["Glyma.TARGET"], _pretty(e_eval)
    assert e_eval["missing_amplicon_gene_mapping"]["target"] == [
        {
            "subject_id": "XM_UNKNOWN",
            "gene_id": "XM_UNKNOWN",
            "glyma_id": "",
        }
    ], _pretty(e_eval)

    assert e_eval["target_expr_by_tissue"] == {"root": 3.0}, _pretty(e_eval)


def test_ratio_is_zero_when_target_missing():
    pair = {
        "target_amplicons": [],
        "off_targets": {
            "near": [{"subject_id": "XM_NEAR", "gene_id": "XM_NEAR", "glyma_id": "", "length": 720}],
            "far": [],
        },
    }

    client = FakeExpressionClient(
        {
            "Glyma.NEAR": {
                "raw_expr": {"ROOT": 2.0},
                "raw_to_canonical": {"ROOT": "root"},
                "canonical_expr": {"root": 2.0},
            }
        }
    )

    e_eval = build_expression_summary_for_pair(
        pair,
        expression_client=client,
        resolver=fake_resolver,
        use_cache=False,
    )

    assert e_eval["target_expr_max"] == 0.0, _pretty(e_eval)
    assert e_eval["target_expr_sum"] == 0.0, _pretty(e_eval)
    assert e_eval["target_vs_near_ratio"] == 0.0, _pretty(e_eval)
    assert e_eval["target_vs_far_ratio"] == 0.0, _pretty(e_eval)


def test_evaluate_blast_results_attaches_e_eval():
    blast_results = [make_fake_pair()]

    client = FakeExpressionClient(
        {
            "Glyma.TARGET": {
                "raw_expr": {"84HAS RH": 5.0},
                "raw_to_canonical": {"84HAS RH": "root_hair"},
                "canonical_expr": {"root_hair": 5.0},
            },
            "Glyma.NEAR": {
                "raw_expr": {"ROOT": 2.1},
                "raw_to_canonical": {"ROOT": "root"},
                "canonical_expr": {"root": 2.1},
            },
            "Glyma.FAR": {
                "raw_expr": {"FLOWER": 0.8},
                "raw_to_canonical": {"FLOWER": "flower"},
                "canonical_expr": {"flower": 0.8},
            },
        }
    )

    out = evaluate_blast_results(
        blast_results,
        expression_client=client,
        resolver=fake_resolver,
        use_cache=False,
    )

    assert len(out) == 1
    assert "e_eval" in out[0], _pretty(out[0])
    assert out[0]["e_eval"]["target_expr_by_tissue"] == {"root_hair": 5.0}, _pretty(out[0])