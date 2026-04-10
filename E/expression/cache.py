from __future__ import annotations

import json
import sqlite3
import time
from typing import Any, Dict, Optional


class ExpressionCache:
    """
    SQLite cache for expression matrices.

    Key: (source, gene_id)
    Value: JSON string of expr_by_condition + metadata
    """
    def __init__(self, db_path: str = "expression_cache.sqlite"):
        self.db_path = db_path
        self._init_db()

    def _init_db(self) -> None:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute("""
        CREATE TABLE IF NOT EXISTS expr_cache (
            source TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            payload_json TEXT NOT NULL,
            timestamp REAL NOT NULL,
            PRIMARY KEY (source, gene_id)
        )
        """)
        conn.commit()
        conn.close()

    def get(self, source: str, gene_id: str, ttl_seconds: int = 7 * 24 * 3600) -> Optional[Dict[str, Any]]:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute(
            "SELECT payload_json, timestamp FROM expr_cache WHERE source=? AND gene_id=?",
            (source, gene_id),
        )
        row = cur.fetchone()
        conn.close()
        if not row:
            return None
        payload_json, ts = row
        if ttl_seconds > 0 and (time.time() - float(ts)) > ttl_seconds:
            return None
        try:
            return json.loads(payload_json)
        except Exception:
            return None

    def set(self, source: str, gene_id: str, payload: Dict[str, Any]) -> None:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute(
            "INSERT OR REPLACE INTO expr_cache (source, gene_id, payload_json, timestamp) VALUES (?, ?, ?, ?)",
            (source, gene_id, json.dumps(payload, ensure_ascii=False), time.time()),
        )
        conn.commit()
        conn.close()
