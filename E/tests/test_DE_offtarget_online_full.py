# -*- coding: utf-8 -*-
"""
test_DE_offtarget_online_full.py

用途：
- 必须联网：D 用 NCBI BLAST（Bio.Blast.NCBIWWW），E 用 resolver + expression
- 全量跑完 C 产生的 pairs（你的例子里是 9 或 23，以实际 build_pairs 输出为准）
- 保持并发策略（ThreadPoolExecutor）
- D 阶段加入细粒度进度输出：
  - 每个 primer 完成数 / 耗时 / hits 数 / 失败数
- D 完成后自动做“诊断输出”（无需你手动加任何东西）：
  - 对前 3 个 pair，打印 forward/reverse 的 subject_id 数量、交集大小
  - 打印交集的示例 subject_id（用于判断为什么 amplicons=0）
  - 如果交集>0 但 amplicons=0，提示可能是 HSP 选择/坐标/长度过滤导致

运行：
在项目根目录（C:\\Users\\42006\\Desktop\\new）：
  python test_DE_offtarget_online_full.py
"""

from __future__ import annotations

import time
import threading
from typing import Dict, List, Any, Set, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

from scan_primers import scan_primers
from filter_primers import filter_primers
from pair_builder import build_pairs

from blast_client_d import BlastClientD, BlastHit
from E.evaluate import evaluate_d_output_e


# 你给的新序列
SEQ = (
"ATGGAAGGATCATCAGGTGTGAGGAAAGGCGCATGGAGTCAATTTGAAGATGACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTCAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTAAGATGGTTGAACTATCTAAAACCAAATATCAAGCGGGGAGATTTTAATGAAGATGAGGTTGATTTGATGATCAGATTGCATAAGCTTTTGGGAAACAGATGGTCCCTAATTGCTGGGAGACTTCCGGGAAGAACTTCAAATGATGTTAAGAACTACTGGAACGCCTACATGCGCCGTAAAGTTCACTCTCACAACAAAGACAACAAAATAAAGAAGCAAGAGACCAAATCAACAGTGAAACCCCACGAAGTGATTAAGCCTATTCCTCGAGTTTTAACAAAAACATCCCCATGGTTGCAACGGAAATTCATTAATAGTCCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACATCATCAGAGAATTGGTGGGAGACATTGTTAGCTGACAAGGAAGATAATGCAGTATTTAACAACAACAACAACACGTGCTTCTTTGGTGGGGTGCATGGAGAGCTTAACCTTTGGAACGAAGAGCTTACTTCAATTGATTTTGACTTTGTTACACAAGGTGAACTTGGAGGGATTGTGAACTTTCAAAATCCATAG"
)


class BlastClientDWithProgress(BlastClientD):
    """
    只覆写 unique primer 的 batch BLAST 这一步，加进度输出；
    其它逻辑（配对、分桶、输出 schema）完全复用原实现。
    """

    def _batch_blast_unique_primers(self, primer_seqs: Set[str], max_workers: int) -> Dict[str, List[BlastHit]]:
        seq_to_hits: Dict[str, List[BlastHit]] = {}

        if self.dummy:
            for s in primer_seqs:
                seq_to_hits[s] = self._blast_short(s)
            return seq_to_hits

        total = len(primer_seqs)
        done = 0
        fail = 0
        lock = threading.Lock()
        t0 = time.perf_counter()

        def short(s: str) -> str:
            if not s:
                return ""
            return s[:10] + "..." + s[-6:] if len(s) > 20 else s

        def task(seq: str) -> List[BlastHit]:
            start = time.perf_counter()
            ok = True
            err = None
            try:
                hits = self._blast_short(seq)
            except Exception as e:
                hits = []
                ok = False
                err = repr(e)

            dt = time.perf_counter() - start
            nonlocal done, fail
            with lock:
                done += 1
                if not ok:
                    fail += 1
                elapsed = time.perf_counter() - t0
                status = "OK" if ok else "FAIL"
                print(
                    f"[D] {done}/{total} {status} "
                    f"primer={short(seq)} "
                    f"time={dt:.1f}s hits={len(hits)} "
                    f"fail={fail} elapsed={elapsed:.1f}s"
                )
                if err and not ok:
                    print(f"[D]   error={err[:200]}")
            return hits

        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            future_to_seq = {ex.submit(task, s): s for s in primer_seqs}
            for fut in as_completed(future_to_seq):
                seq = future_to_seq[fut]
                try:
                    seq_to_hits[seq] = fut.result()
                except Exception:
                    seq_to_hits[seq] = []

        elapsed = time.perf_counter() - t0
        print(f"[D] batch finished: total={total} done={done} fail={fail} elapsed={elapsed:.1f}s")
        return seq_to_hits


def _subjects_from_hits(hit_dicts: List[Dict[str, Any]]) -> Set[str]:
    out: Set[str] = set()
    for h in hit_dicts or []:
        sid = (h.get("subject_id") or "").strip()
        if sid:
            out.add(sid)
    return out


def _print_pair_subject_diagnostics(d_out: List[Dict[str, Any]], n_pairs: int = 3) -> None:
    """
    自动诊断：为什么 amplicons=0。
    - 看 forward/reverse 命中 subject_id 是否有交集
    - 如果交集为 0：说明 subject_id 对不上（最常见）
    - 如果交集 > 0 但 amplicons=0：说明坐标/长度/取 HSP 的策略在过滤
    """
    print("\n[DIAG] subject_id intersection diagnostics (first pairs)")

    for i, pr in enumerate(d_out[:n_pairs]):
        hits = pr.get("hits") or {}
        f_hits = hits.get("forward") or []
        r_hits = hits.get("reverse") or []

        f_subjects = _subjects_from_hits(f_hits)
        r_subjects = _subjects_from_hits(r_hits)
        inter = f_subjects & r_subjects

        amps = pr.get("amplicons") or []
        print(f"\n[DIAG] Pair {pr.get('pair_index', i)}")
        print(f"[DIAG]  forward_hits={len(f_hits)} reverse_hits={len(r_hits)} amplicons={len(amps)}")
        print(f"[DIAG]  f_subjects={len(f_subjects)} r_subjects={len(r_subjects)} intersection={len(inter)}")

        if len(inter) == 0:
            # 打印少量示例，帮助你确认 subject_id 形态
            ex_f = list(sorted(f_subjects))[:5]
            ex_r = list(sorted(r_subjects))[:5]
            print("[DIAG]  intersection=0 => forward/reverse 的 subject_id 没对上（最常见原因）")
            print("[DIAG]  example forward subject_id:", ex_f)
            print("[DIAG]  example reverse subject_id:", ex_r)
        else:
            ex_i = list(sorted(inter))[:10]
            print("[DIAG]  example intersection subject_id:", ex_i)
            if len(amps) == 0:
                print("[DIAG]  intersection>0 但 amplicons=0 => 更可能是坐标/长度/取 HSP(只取hsps[0])导致过滤")


def main() -> None:
    sequence = SEQ.replace("\n", "").replace(" ", "").upper()

    # ---------------- A ----------------
    ok, primers, info = scan_primers(sequence, size_range=(18, 25))
    assert ok, f"A failed: {info}"

    # ---------------- B ----------------
    fwd = filter_primers(
        primers["forward"],
        min_tm=50.0,
        max_tm=65.0,
        max_hairpin_dG=-3.0,
        max_dimer_dG=-9.0,
    )
    rev = filter_primers(
        primers["reverse"],
        min_tm=50.0,
        max_tm=65.0,
        max_hairpin_dG=-3.0,
        max_dimer_dG=-9.0,
    )
    assert len(fwd) > 0 and len(rev) > 0, "B produced empty primer lists"

    # ---------------- C ----------------
    tl = len(sequence)
    pairs = build_pairs(
        fwd,
        rev,
        target_length=tl,
        max_tm_diff=2.0,
        max_heterodimer_dG=0.0,
    )
    assert len(pairs) > 0, "C produced no primer pairs"

    print(f"Total pairs generated: {len(pairs)}")
    print(f"[ABC] fwd={len(fwd)} rev={len(rev)} pairs={len(pairs)} target_length={tl}")

    # ---------------- D (online, async concurrency + progress) ----------------
    # 仍然是异步并发；为了更容易产生候选命中，hitlist_size 设置为 50。
    # 你可以把 hitlist_size 改到 100，但 NCBI 压力会更大。
    max_workers = 2

    client = BlastClientDWithProgress(
        organism="Glycine max",
        db="refseq_rna",
        hitlist_size=50,
        expect=1000.0,
        word_size=7,
        megablast=False,
        max_retries=3,
        backoff_base=2.0,
        sleep_between_requests=0.4,
        dummy=False,
    )

    t0 = time.perf_counter()
    d_out = client.process_pairs(
        primer_pairs=pairs,
        max_workers=max_workers,
        target_length=tl,
        near_bp=200,
        target_subject_id=None,
        min_amplicon_len=50,
        max_amplicon_len=3000,
        require_strict_orientation=False,
        min_identity_min=None,
    )
    dt = time.perf_counter() - t0
    print(f"[D] pairs processed: {len(d_out)} elapsed={dt:.1f}s")

    # 汇总 D 的 amplicons / off-target 情况
    near_total = 0
    far_total = 0
    amp_total = 0
    for pr in d_out:
        amp_total += len(pr.get("amplicons") or [])
        off = pr.get("off_targets") or {}
        near_total += len(off.get("near") or [])
        far_total += len(off.get("far") or [])
    print(f"[D] total_amplicons={amp_total} off_near={near_total} off_far={far_total}")

    # 自动诊断：为什么 total_amplicons=0
    _print_pair_subject_diagnostics(d_out, n_pairs=3)

    # ---------------- E (online) ----------------
    print("[E] entering evaluate...")
    t1 = time.perf_counter()
    de_out = evaluate_d_output_e(
        d_out,
        expression_cache_db_path="expression_cache.sqlite",
        resolver_cache_db_path="resolver_cache.sqlite",
        use_cache=True,
        enable_ensembl_fallback=True,
        conekt_sleep_s=0.2,
        conekt_timeout=25,
    )
    dt_e = time.perf_counter() - t1
    print(f"[E] evaluate finished: results={len(de_out)} elapsed={dt_e:.1f}s")

    # ---------------- 打印前 3 个 pair 的关键字段 ----------------
    for pr in de_out[:3]:
        ee = pr.get("e_eval") or {}
        off_ids = ee.get("off_targets_glyma_ids") or {}
        print("\n--- Pair", pr.get("pair_index"), "---")
        print("target_subject_id:", pr.get("target_subject_id"))
        print("target_glyma_id:", ee.get("target_glyma_id"))
        print("off_near_ids:", len(off_ids.get("near") or []), "off_far_ids:", len(off_ids.get("far") or []))
        print("expr:", ee.get("expr"))
        print("missing:", ee.get("missing"))

    if amp_total == 0:
        print("\n[RESULT] amplicons=0 => D 配对阶段未形成扩增产物。")
        print("         重点看上面的 [DIAG]：")
        print("         - 如果 intersection=0：应优先修复 subject_id 归一化（让 F/R 落到同一个 subject key）")
        print("         - 如果 intersection>0 但 amplicons=0：应修复 HSP 选择/坐标与长度过滤逻辑")
    else:
        print("\n[RESULT] amplicons>0 => D 配对已覆盖；E 的 off-target 聚合将被真实覆盖。")


if __name__ == "__main__":
    main()
