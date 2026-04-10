from E.expression.client import ExpressionClient

client = ExpressionClient()

gene = "Glyma.09G235300"

ge = client.fetch_gene(gene)

print("\n=== EXPRESSION TEST ===")

print("gene:", gene)
print("ok:", ge.ok)
print("source:", ge.source)
print("expr:", ge.raw_expr)
print("notes:", ge.notes)