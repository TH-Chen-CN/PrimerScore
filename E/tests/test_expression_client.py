from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from E.expression.client import ExpressionClient  # noqa: E402


def _pretty(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True, default=str)


class FakeCache:
    def __init__(self):
        self.store = {}

    def get(self, source, prefixed, ttl_seconds=None):
        return self.store.get((source, prefixed))

    def set(self, source, prefixed, payload):
        self.store[(source, prefixed)] = payload


class FakeConektClient:
    def __init__(self, payload_by_gene):
        self.payload_by_gene = payload_by_gene

    def fetch_gene_expression(self, prefixed_gene):
        payload = self.payload_by_gene.get(prefixed_gene)
        if payload is None:
            return SimpleNamespace(
                ok=False,
                expr_by_condition={},
                notes="conekt miss",
                source_url="https://conekt.example/miss",
            )
        return SimpleNamespace(
            ok=True,
            expr_by_condition=payload,
            notes="conekt ok",
            source_url="https://conekt.example/hit",
        )


class FakeSoyOmicsClient:
    def __init__(self, payload_by_gene):
        self.payload_by_gene = payload_by_gene

    def fetch_gene_expression(self, glyma_gene):
        payload = self.payload_by_gene.get(glyma_gene)
        if payload is None:
            return SimpleNamespace(
                ok=False,
                expr_by_condition={},
                notes="soyomics miss",
                source_url="https://soyomics.example/miss",
            )
        return SimpleNamespace(
            ok=True,
            expr_by_condition=payload,
            notes="soyomics ok",
            source_url="https://soyomics.example/hit",
        )


def make_client(
    conekt_payload=None,
    soyomics_payload=None,
    preferred="integrated",
    fallback=True,
):
    client = ExpressionClient(
        cache_db_path=":memory:",
        preferred_source=preferred,
        enable_soyomics_fallback=fallback,
    )
    client.cache = FakeCache()
    client.conekt = FakeConektClient(conekt_payload or {})
    client.soyomics = FakeSoyOmicsClient(soyomics_payload or {})
    return client


def test_fetch_gene_prefers_conekt():
    client = make_client(
        conekt_payload={
            "glyma.Wm82.gnm2.ann1.Glyma.19G020600": {
                "84HAS RH": 5.0,
                "ROOT": 1.0,
            }
        },
        soyomics_payload={
            "Glyma.19G020600": {
                "stem": 9.0,
            }
        },
        preferred="conekt_lis",
        fallback=True,
    )

    ge = client.fetch_gene("Glyma.19G020600", use_cache=False)

    assert ge.ok is True, _pretty(ge.__dict__)
    assert ge.source == "conekt_lis", _pretty(ge.__dict__)
    assert ge.canonical_expr["root_hair"] == 5.0, _pretty(ge.__dict__)
    assert ge.canonical_expr["root"] == 1.0, _pretty(ge.__dict__)


def test_fetch_gene_fallback_to_soyomics():
    client = make_client(
        conekt_payload={},
        soyomics_payload={
            "Glyma.09G235300": {
                "stem": 8.5,
                "root": 0.0,
                "flower": 0.0,
            }
        },
        preferred="conekt_lis",
        fallback=True,
    )

    ge = client.fetch_gene("Glyma.09G235300", use_cache=False)

    assert ge.ok is True, _pretty(ge.__dict__)
    assert ge.source == "soyomics", _pretty(ge.__dict__)
    assert ge.canonical_expr["hypocotyl"] == 8.5, _pretty(ge.__dict__)


def test_fetch_gene_prefers_soyomics():
    client = make_client(
        conekt_payload={
            "glyma.Wm82.gnm2.ann1.Glyma.09G235300": {
                "24HA1_IN_RH": 0.83,
            }
        },
        soyomics_payload={
            "Glyma.09G235300": {
                "stem": 12.0,
                "root": 0.0,
            }
        },
        preferred="soyomics",
        fallback=True,
    )

    ge = client.fetch_gene("Glyma.09G235300", use_cache=False)

    assert ge.ok is True, _pretty(ge.__dict__)
    assert ge.source == "soyomics", _pretty(ge.__dict__)
    assert ge.canonical_expr["hypocotyl"] == 12.0, _pretty(ge.__dict__)


def test_fetch_gene_no_fallback():
    client = make_client(
        conekt_payload={},
        soyomics_payload={
            "Glyma.09G235300": {"stem": 5.0}
        },
        preferred="conekt_lis",
        fallback=False,
    )

    ge = client.fetch_gene("Glyma.09G235300", use_cache=False)

    assert ge.ok is False, _pretty(ge.__dict__)
    assert ge.raw_expr == {}, _pretty(ge.__dict__)


def test_cache_roundtrip():
    client = make_client(
        conekt_payload={
            "glyma.Wm82.gnm2.ann1.Glyma.19G020600": {
                "84HAS RH": 5.0,
            }
        },
        preferred="conekt_lis",
        fallback=True,
    )

    ge1 = client.fetch_gene("Glyma.19G020600", use_cache=True)
    ge2 = client.fetch_gene("Glyma.19G020600", use_cache=True)

    assert ge1.raw_expr == ge2.raw_expr
    assert ge1.canonical_expr == ge2.canonical_expr
    assert ge1.source == ge2.source


def test_integrated_source_merges_soyomics_and_conekt():
    client = make_client(
        conekt_payload={
            "glyma.Wm82.gnm2.ann1.Glyma.09G235300": {
                "24HA1_IN_RH": 0.83,
                "ROOT": 0.1,
            }
        },
        soyomics_payload={
            "Glyma.09G235300": {
                "stem": 1.6,
                "flower": 0.0,
                "root": 0.0,
            }
        },
        preferred="integrated",
        fallback=True,
    )

    ge = client.fetch_gene("Glyma.09G235300", use_cache=False)

    assert ge.ok is True, _pretty(ge.__dict__)
    assert ge.source == "integrated", _pretty(ge.__dict__)
    assert ge.canonical_expr["hypocotyl"] == 1.6, _pretty(ge.__dict__)
    assert ge.canonical_expr["root_hair"] == 0.83, _pretty(ge.__dict__)
    assert "soyomics" in ge.by_source, _pretty(ge.__dict__)
    assert "conekt_lis" in ge.by_source, _pretty(ge.__dict__)