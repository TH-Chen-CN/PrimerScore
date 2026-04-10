"""
CoNekT expression client (SoyBase/LIS).

Goal:
- Input: prefixed gene id, e.g. "glyma.Wm82.gnm2.ann1.Glyma.19G020600"
- Output: dict[raw_condition] = mean_expression(float)

How it works:
1) Open gene page: /sequence/find_forLIS/<gene_id>
2) Extract the "Download" link under "Expression Profile"
3) Download plain text from /profile/download/plot/<profile_id>
4) Parse: columns = condition mean min max

Notes:
- This is perfect for E because expression values are internally comparable across conditions for the same dataset.
"""

from __future__ import annotations

import re
import time
import requests
from dataclasses import dataclass
from typing import Dict, Optional, Tuple


DEFAULT_TIMEOUT = 20
DEFAULT_UA = "primer-e-expression/0.1 (+https://example.local)"
CONEKT_BASE = "https://conekt.legumeinfo.org"


@dataclass
class FetchResult:
    ok: bool
    expr_by_condition: Dict[str, float]
    notes: str
    source_url: str


class ConektExpressionClient:
    def __init__(
        self,
        timeout: int = DEFAULT_TIMEOUT,
        user_agent: str = DEFAULT_UA,
        sleep_s: float = 0.0,
    ):
        self.timeout = timeout
        self.headers = {"User-Agent": user_agent}
        self.sleep_s = sleep_s

    def fetch_gene_expression(self, prefixed_gene_id: str) -> FetchResult:
        """
        Return expression means by condition for a single gene.
        """
        gene_url = f"{CONEKT_BASE}/sequence/find_forLIS/{prefixed_gene_id}"

        try:
            if self.sleep_s:
                time.sleep(self.sleep_s)

            r = requests.get(gene_url, headers=self.headers, timeout=self.timeout)
            if r.status_code != 200:
                return FetchResult(
                    ok=False,
                    expr_by_condition={},
                    notes=f"gene page HTTP {r.status_code}",
                    source_url=gene_url,
                )

            html = r.text

            # Find the download link for expression profile
            # Example: https://conekt.legumeinfo.org/profile/download/plot/27201
            m = re.search(r'href="(/profile/download/plot/\d+)"', html)
            if not m:
                return FetchResult(
                    ok=False,
                    expr_by_condition={},
                    notes="expression download link not found in gene page",
                    source_url=gene_url,
                )

            download_path = m.group(1)
            download_url = f"{CONEKT_BASE}{download_path}"

            if self.sleep_s:
                time.sleep(self.sleep_s)

            t = requests.get(download_url, headers=self.headers, timeout=self.timeout)
            if t.status_code != 200:
                return FetchResult(
                    ok=False,
                    expr_by_condition={},
                    notes=f"download HTTP {t.status_code}",
                    source_url=download_url,
                )

            expr = self._parse_download_plot_text(t.text)
            if not expr:
                return FetchResult(
                    ok=False,
                    expr_by_condition={},
                    notes="download parsed empty",
                    source_url=download_url,
                )

            return FetchResult(
                ok=True,
                expr_by_condition=expr,
                notes="ok",
                source_url=download_url,
            )

        except requests.RequestException as e:
            return FetchResult(
                ok=False,
                expr_by_condition={},
                notes=f"requests error: {e.__class__.__name__}: {e}",
                source_url=gene_url,
            )
        except Exception as e:
            return FetchResult(
                ok=False,
                expr_by_condition={},
                notes=f"unexpected error: {e.__class__.__name__}: {e}",
                source_url=gene_url,
            )

    @staticmethod
    def _parse_download_plot_text(text: str) -> Dict[str, float]:
        """
        Input looks like one long line (sometimes):
        condition mean min max 12HA1_IN_RH 0.0 0.0 0.0 Stacey_Apical_Meristem 1.1 1.1 1.1 ...

        We parse tokens in groups of 4 after the header.
        """
        if not text:
            return {}

        tokens = text.strip().split()
        if len(tokens) < 8:
            return {}

        # Must start with: condition mean min max
        if tokens[0:4] != ["condition", "mean", "min", "max"]:
            # Sometimes may have newlines; try normalize
            pass

        # Find header index safely
        try:
            idx = tokens.index("condition")
            if tokens[idx:idx + 4] != ["condition", "mean", "min", "max"]:
                return {}
            start = idx + 4
        except ValueError:
            return {}

        out: Dict[str, float] = {}
        i = start
        while i + 3 < len(tokens):
            cond = tokens[i]
            mean_s = tokens[i + 1]
            # tokens[i+2] min, tokens[i+3] max
            try:
                mean_v = float(mean_s)
            except ValueError:
                # If parse breaks, stop to avoid garbage.
                break
            out[cond] = mean_v
            i += 4

        return out


def ensure_prefixed_glyma(glyma_id: str, prefix: str = "glyma.Wm82.gnm2.ann1.") -> str:
    """
    glyma_id: "Glyma.19G020600" -> "glyma.Wm82.gnm2.ann1.Glyma.19G020600"
    If already prefixed, keep.
    """
    if glyma_id.startswith("glyma."):
        return glyma_id
    return f"{prefix}{glyma_id}"
