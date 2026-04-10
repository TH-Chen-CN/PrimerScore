from __future__ import annotations

import json
from pathlib import Path
import sys

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from E.expression.canonical import (  # noqa: E402
    CANONICAL_TISSUES,
    canonicalize_condition,
    canonical_condition,
)


def _pretty(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True)


def test_vocab_contains_soyomics_expansion():
    expected_subset = {
        "cotyledon",
        "hypocotyl",
        "shoot_apex",
        "seed",
        "axillary_bud",
    }
    assert expected_subset.issubset(CANONICAL_TISSUES), (
        f"Missing soyomics tissues.\nActual: {sorted(CANONICAL_TISSUES)}"
    )


def test_old_atlas_root_hair_mapping_still_works():
    raw = "12HA1_IN_RH"
    result = canonicalize_condition(raw, source="conekt_lis")
    assert result["canonical_tissue"] == "root_hair", _pretty(result)


def test_soyomics_cotyledon_mapping():
    raw = "cotyledon-1"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "cotyledon", _pretty(result)


def test_soyomics_stem_maps_to_hypocotyl():
    raw = "stem-1"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "hypocotyl", _pretty(result)


def test_soyomics_leafbud_maps_to_shoot_apex():
    raw = "leafbud-2"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "shoot_apex", _pretty(result)


def test_soyomics_seed_mapping():
    raw = "seed-3"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "seed", _pretty(result)


def test_soyomics_axillary_bud_mapping():
    raw = "axillary_bud"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "axillary_bud", _pretty(result)


def test_soyomics_axillary_bud_space_variant():
    raw = "axillary bud"
    result = canonicalize_condition(raw, source="soyomics")
    assert result["canonical_tissue"] == "axillary_bud", _pretty(result)


def test_canonical_condition_helper():
    assert canonical_condition("stem-1", source="soyomics") == "hypocotyl"
    assert canonical_condition("12HA1_IN_RH", source="conekt_lis") == "root_hair"