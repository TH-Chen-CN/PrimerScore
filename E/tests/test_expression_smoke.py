from __future__ import annotations

from E.expression.client import ExpressionClient


def main():
    gene = "Glyma.19G020600"
    client = ExpressionClient(
        allowed_canon_conditions={"stem", "leaf", "root", "seed", "flower", "pod", "nodule"},
        conekt_sleep_s=0.2,
    )

    ge = client.fetch_gene(gene, use_cache=True)

    print("Gene:", ge.glyma_gene)
    print("Prefixed:", ge.prefixed_gene)
    print("OK:", ge.ok)
    print("Notes:", ge.notes)
    print("Source URL:", ge.source_url)
    print("Raw n_conditions:", len(ge.raw_expr))
    print("Canonical n_conditions:", len(ge.canonical_expr))

    print("\n[Raw conditions -> canonical]")
    for k in list(ge.raw_expr.keys())[:20]:
        print(f"- {k} -> {ge.raw_to_canonical.get(k)}  value={ge.raw_expr[k]}")

    print("\n[Canonical expr (aggregated by max)]")
    for k, v in sorted(ge.canonical_expr.items()):
        print(f"- {k}: {v}")


if __name__ == "__main__":
    main()
