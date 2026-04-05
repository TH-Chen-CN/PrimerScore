from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any, Set
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import time

from Bio.Blast import NCBIWWW, NCBIXML


# ----------------------------
# Data models (stable schema)
# ----------------------------

@dataclass(frozen=True)
class BlastHit:
    subject_id: str
    title: str
    identity: float
    e_value: float
    align_len: int
    sbjct_start: int
    sbjct_end: int
    strand: str  # "+" or "-"


@dataclass(frozen=True)
class Amplicon:
    subject_id: str
    gene_id: str       # stable-ish label (prefer subject_id)
    glyma_id: str      # extracted Glyma.* if available, else ""
    length: int
    fwd: BlastHit
    rev: BlastHit
    identity_min: float
    e_value_max: float


@dataclass
class PairResult:
    pair_index: int
    forward: str
    reverse: str
    product_size: Optional[int]
    tm_diff: Optional[float]
    heterodimer_dG: Optional[float]

    target_subject_id: Optional[str]  # NEW helper

    hits: Dict[str, List[Dict[str, Any]]]
    amplicons: List[Dict[str, Any]]

    target_amplicons: List[Dict[str, Any]]
    off_targets: Dict[str, List[Dict[str, Any]]]

    status: str
    error: Optional[str] = None


# ----------------------------
# D module implementation
# ----------------------------

class BlastClientD:
    def __init__(
        self,
        organism: str,
        db: str = "refseq_rna",
        hitlist_size: int = 10,
        expect: float = 1000.0,
        word_size: int = 7,
        megablast: bool = False,
        max_retries: int = 3,
        backoff_base: float = 2.0,
        sleep_between_requests: float = 0.0,
        dummy: bool = False,
    ):
        self.organism = organism
        self.db = db
        self.hitlist_size = hitlist_size
        self.expect = expect
        self.word_size = word_size
        self.megablast = megablast
        self.max_retries = max_retries
        self.backoff_base = backoff_base
        self.sleep_between_requests = sleep_between_requests
        self.dummy = dummy

        self._cache: Dict[Tuple[str, str, str], List[BlastHit]] = {}
        self._cache_lock = threading.Lock()

    # --------- Public API ---------

    def process_pairs(
        self,
        primer_pairs: List[Dict[str, Any]],
        max_workers: int = 6,
        target_length: Optional[int] = None,
        near_bp: int = 30,
        target_subject_id: Optional[str] = None,
        min_amplicon_len: int = 50,
        max_amplicon_len: int = 3000,
        require_strict_orientation: bool = True,
        # NEW: filter amplicons by identity_min (amplicon-level)
        min_identity_min: Optional[float] = None,
    ) -> List[Dict[str, Any]]:
        if max_workers < 1:
            max_workers = 1

        unique_primers = self._collect_unique_primers(primer_pairs)
        seq_to_hits = self._batch_blast_unique_primers(unique_primers, max_workers=max_workers)

        results: List[PairResult] = []
        for idx, pair in enumerate(primer_pairs):
            results.append(
                self._build_pair_result(
                    idx=idx,
                    pair=pair,
                    seq_to_hits=seq_to_hits,
                    target_length=target_length,
                    near_bp=near_bp,
                    target_subject_id=target_subject_id,
                    min_amplicon_len=min_amplicon_len,
                    max_amplicon_len=max_amplicon_len,
                    require_strict_orientation=require_strict_orientation,
                    min_identity_min=min_identity_min,  # NEW
                )
            )

        return [asdict(r) for r in results]

    # --------- Step A: scatter/gather BLAST ---------

    def _collect_unique_primers(self, primer_pairs: List[Dict[str, Any]]) -> Set[str]:
        uniq: Set[str] = set()
        for pair in primer_pairs:
            f = pair.get("forward", {}).get("sequence")
            r = pair.get("reverse", {}).get("sequence")
            if isinstance(f, str) and f:
                uniq.add(f)
            if isinstance(r, str) and r:
                uniq.add(r)
        return uniq

    def _batch_blast_unique_primers(self, primer_seqs: Set[str], max_workers: int) -> Dict[str, List[BlastHit]]:
        seq_to_hits: Dict[str, List[BlastHit]] = {}

        if self.dummy:
            for s in primer_seqs:
                seq_to_hits[s] = self._blast_short(s)
            return seq_to_hits

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            future_to_seq = {ex.submit(self._blast_short, s): s for s in primer_seqs}
            for fut in as_completed(future_to_seq):
                seq = future_to_seq[fut]
                try:
                    seq_to_hits[seq] = fut.result()
                except Exception:
                    seq_to_hits[seq] = []

        return seq_to_hits

    # --------- BLAST core ---------

    @staticmethod
    def _hsp_subject_strand(hsp) -> str:
        frame = getattr(hsp, "frame", None)
        if isinstance(frame, (tuple, list)) and len(frame) >= 2:
            sbj = frame[1]
            if isinstance(sbj, int):
                return "+" if sbj >= 0 else "-"

        strand = getattr(hsp, "strand", None)
        if isinstance(strand, (tuple, list)) and len(strand) >= 2:
            sbj = strand[1]
            if isinstance(sbj, str):
                s = sbj.lower()
                if "minus" in s:
                    return "-"
                if "plus" in s:
                    return "+"

        return "+"

    def _blast_short(self, seq: str) -> List[BlastHit]:
        cache_key = (seq, self.organism, self.db)
        with self._cache_lock:
            if cache_key in self._cache:
                return self._cache[cache_key]

        if self.dummy:
            hits = self._dummy_hits(seq)
            with self._cache_lock:
                self._cache[cache_key] = hits
            return hits

        last_err: Optional[Exception] = None
        for attempt in range(self.max_retries):
            try:
                if self.sleep_between_requests > 0:
                    time.sleep(self.sleep_between_requests)

                handle = NCBIWWW.qblast(
                    program="blastn",
                    database=self.db,
                    sequence=seq,
                    url_base="https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                    megablast=self.megablast,
                    expect=self.expect,
                    word_size=self.word_size,
                    hitlist_size=self.hitlist_size,
                    service="plain",
                    entrez_query=f"{self.organism}[organism]",
                    format_type="XML",
                )

                record = NCBIXML.read(handle)
                hits: List[BlastHit] = []

                for aln in getattr(record, "alignments", []):
                    hsps = getattr(aln, "hsps", None)
                    if not hsps:
                        continue
                    hsp = hsps[0]

                    subject_id = getattr(aln, "accession", "") or self._extract_subject_id(getattr(aln, "title", ""))
                    strand = self._hsp_subject_strand(hsp)

                    hits.append(BlastHit(
                        subject_id=subject_id,
                        title=(aln.title[:200] if getattr(aln, "title", None) else ""),
                        identity=round(hsp.identities * 100.0 / max(1, hsp.align_length), 2),
                        e_value=float(hsp.expect),
                        align_len=int(hsp.align_length),
                        sbjct_start=int(hsp.sbjct_start),
                        sbjct_end=int(hsp.sbjct_end),
                        strand=strand,
                    ))

                with self._cache_lock:
                    self._cache[cache_key] = hits
                return hits

            except Exception as e:
                last_err = e
                time.sleep(self.backoff_base ** attempt)

        with self._cache_lock:
            self._cache[cache_key] = []
        if last_err:
            raise last_err
        return []

    # --------- Step B/C: pairing + bucketing ---------

    def _build_pair_result(
        self,
        idx: int,
        pair: Dict[str, Any],
        seq_to_hits: Dict[str, List[BlastHit]],
        target_length: Optional[int],
        near_bp: int,
        target_subject_id: Optional[str],
        min_amplicon_len: int,
        max_amplicon_len: int,
        require_strict_orientation: bool,
        # NEW
        min_identity_min: Optional[float],
    ) -> PairResult:
        try:
            fwd_seq = pair["forward"]["sequence"]
            rev_seq = pair["reverse"]["sequence"]

            product_size = pair.get("product_size", None)
            tm_diff = pair.get("tm_diff", None)
            heterodimer_dG = pair.get("heterodimer_dG", None)

            tl = target_length
            if tl is None and isinstance(product_size, int):
                tl = product_size

            fwd_hits = seq_to_hits.get(fwd_seq, [])
            rev_hits = seq_to_hits.get(rev_seq, [])

            amplicons = self._pair_hits_into_amplicons_strict(
                fwd_hits=fwd_hits,
                rev_hits=rev_hits,
                min_len=min_amplicon_len,
                max_len=max_amplicon_len,
                require_strict_orientation=require_strict_orientation,
            )

            # NEW: filter amplicons by identity_min if requested
            if min_identity_min is not None:
                amplicons = [a for a in amplicons if a.identity_min >= float(min_identity_min)]

            chosen_target_subject = target_subject_id or self._infer_target_subject(amplicons)

            target_amplicons = [a for a in amplicons if a.subject_id == chosen_target_subject] if chosen_target_subject else []
            off_amps = [a for a in amplicons if (chosen_target_subject is None or a.subject_id != chosen_target_subject)]

            near_list: List[Amplicon] = []
            far_list: List[Amplicon] = []
            if tl is not None:
                for a in off_amps:
                    if abs(a.length - tl) <= near_bp:
                        near_list.append(a)
                    else:
                        far_list.append(a)
            else:
                far_list = off_amps

            return PairResult(
                pair_index=idx,
                forward=fwd_seq,
                reverse=rev_seq,
                product_size=product_size,
                tm_diff=tm_diff,
                heterodimer_dG=heterodimer_dG,
                target_subject_id=chosen_target_subject,
                hits={
                    "forward": [asdict(h) for h in fwd_hits],
                    "reverse": [asdict(h) for h in rev_hits],
                },
                amplicons=[asdict(a) for a in amplicons],
                target_amplicons=[asdict(a) for a in target_amplicons],
                off_targets={
                    "near": [asdict(a) for a in near_list],
                    "far": [asdict(a) for a in far_list],
                },
                status="success",
                error=None,
            )

        except Exception as e:
            fwd_seq = pair.get("forward", {}).get("sequence", "")
            rev_seq = pair.get("reverse", {}).get("sequence", "")
            return PairResult(
                pair_index=idx,
                forward=fwd_seq,
                reverse=rev_seq,
                product_size=pair.get("product_size", None),
                tm_diff=pair.get("tm_diff", None),
                heterodimer_dG=pair.get("heterodimer_dG", None),
                target_subject_id=None,
                hits={"forward": [], "reverse": []},
                amplicons=[],
                target_amplicons=[],
                off_targets={"near": [], "far": []},
                status="failed",
                error=str(e),
            )

    def _pair_hits_into_amplicons_strict(
        self,
        fwd_hits: List[BlastHit],
        rev_hits: List[BlastHit],
        min_len: int,
        max_len: int,
        require_strict_orientation: bool = True,
    ) -> List[Amplicon]:
        by_subject_f: Dict[str, List[BlastHit]] = {}
        by_subject_r: Dict[str, List[BlastHit]] = {}

        for h in fwd_hits:
            if h.subject_id:
                by_subject_f.setdefault(h.subject_id, []).append(h)
        for h in rev_hits:
            if h.subject_id:
                by_subject_r.setdefault(h.subject_id, []).append(h)

        amplicons: List[Amplicon] = []

        for subject_id, f_list in by_subject_f.items():
            r_list = by_subject_r.get(subject_id, [])
            if not r_list:
                continue

            for fh in f_list:
                for rh in r_list:
                    if require_strict_orientation:
                        if not (fh.strand == "+" and rh.strand == "-"):
                            continue

                    f_left = min(fh.sbjct_start, fh.sbjct_end)
                    r_right = max(rh.sbjct_start, rh.sbjct_end)

                    if not (f_left < r_right):
                        continue

                    length = (r_right - f_left) + 1
                    if length < min_len or length > max_len:
                        continue

                    glyma_id = self._extract_glyma_id(fh.title) or self._extract_glyma_id(rh.title) or ""
                    gene_id = subject_id or (fh.title[:60].strip() if fh.title else "")

                    amplicons.append(Amplicon(
                        subject_id=subject_id,
                        gene_id=gene_id,
                        glyma_id=glyma_id,
                        length=int(length),
                        fwd=fh,
                        rev=rh,
                        identity_min=min(fh.identity, rh.identity),
                        e_value_max=max(fh.e_value, rh.e_value),
                    ))

        uniq: Dict[Tuple[str, int, int, int], Amplicon] = {}
        for a in amplicons:
            k = (a.subject_id, a.length, a.fwd.sbjct_start, a.rev.sbjct_start)
            if k not in uniq:
                uniq[k] = a
        return list(uniq.values())

    def _infer_target_subject(self, amplicons: List[Amplicon]) -> Optional[str]:
        if not amplicons:
            return None
        ranked = sorted(amplicons, key=lambda a: (-a.identity_min, a.e_value_max, a.length))
        return ranked[0].subject_id

    # --------- Helpers ---------

    def _extract_glyma_id(self, title: str) -> str:
        if not title:
            return ""
        parts = title.split("|")
        for p in parts:
            p = p.strip()
            if p.startswith("Glyma."):
                return p
        return ""

    def _extract_subject_id(self, title: str) -> str:
        if not title:
            return ""
        parts = title.split("|")
        for p in parts:
            p = p.strip()
            if p.startswith(("NM_", "XM_", "NR_", "XR_")):
                return p
        return ""

    def _dummy_hits(self, seq: str) -> List[BlastHit]:
        sid = f"DUMMY_{len(seq)}"
        return [
            BlastHit(
                subject_id=sid,
                title="dummy alignment",
                identity=100.0,
                e_value=0.0,
                align_len=len(seq),
                sbjct_start=100,
                sbjct_end=100 + len(seq) - 1,
                strand="+",
            )
        ]


def process_c_module_output_d(
    c_module_output: List[Dict[str, Any]],
    organism: str,
    db: str = "refseq_rna",
    max_workers: int = 6,
    target_length: Optional[int] = None,
    near_bp: int = 30,
    target_subject_id: Optional[str] = None,
    dummy: bool = False,
    min_amplicon_len: int = 50,
    max_amplicon_len: int = 3000,
    require_strict_orientation: bool = True,
    sleep_between_requests: float = 0.0,
    # NEW
    min_identity_min: Optional[float] = None,
) -> List[Dict[str, Any]]:
    client = BlastClientD(
        organism=organism,
        db=db,
        dummy=dummy,
        sleep_between_requests=sleep_between_requests,
    )
    return client.process_pairs(
        primer_pairs=c_module_output,
        max_workers=max_workers,
        target_length=target_length,
        near_bp=near_bp,
        target_subject_id=target_subject_id,
        min_amplicon_len=min_amplicon_len,
        max_amplicon_len=max_amplicon_len,
        require_strict_orientation=require_strict_orientation,
        min_identity_min=min_identity_min,  # NEW
    )
