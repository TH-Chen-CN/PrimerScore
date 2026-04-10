# E/resolver/mapping_index.py
from __future__ import annotations

import sqlite3
from typing import Optional, Dict


class MappingIndex:
    """
    Lightweight local mapping index:
      transcript_accession -> glyma_id

    This is *not* a test tool; it's a core helper used by resolver.
    """

    def __init__(self, db_path: str = "mapping_index.sqlite"):
        self.db_path = db_path
        self._init_db()

    def _init_db(self) -> None:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute("""
        CREATE TABLE IF NOT EXISTS transcript_to_glyma (
            transcript_id TEXT PRIMARY KEY,
            glyma_id TEXT NOT NULL
        )
        """)
        cur.execute("CREATE INDEX IF NOT EXISTS idx_transcript_id ON transcript_to_glyma(transcript_id)")
        conn.commit()
        conn.close()

    def get(self, transcript_id: str) -> Optional[Dict]:
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        cur.execute(
            "SELECT glyma_id FROM transcript_to_glyma WHERE transcript_id=?",
            (transcript_id,)
        )
        row = cur.fetchone()
        conn.close()
        if not row:
            return None
        return {"glyma_id": row[0]}
