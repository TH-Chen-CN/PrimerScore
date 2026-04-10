# -*- coding: utf-8 -*-
from __future__ import annotations

from scan_primers import scan_primers          # A
from filter_primers import filter_primers      # B
from pair_builder import build_pairs           # C
from blast_client_d import process_c_module_output_d  # D


SEQ = (
"ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTGTTACACAAGGTGAATCTTGGAGTGATTTTCTTCTTGACCTACAAGGCTAG"
)


def assert_pairresult_contract(pr: dict) -> None:
    """
    验证 D 输出的 PairResult 结构契约是否满足 E/F 接入的最低要求。
    """
    required_top_keys = ["amplicons", "target_amplicons", "off_targets", "target_subject_id"]
    for k in required_top_keys:
        assert k in pr, f"PairResult missing key: {k}"

    assert isinstance(pr["amplicons"], list), "amplicons must be a list"
    assert isinstance(pr["target_amplicons"], list), "target_amplicons must be a list"
    assert isinstance(pr["off_targets"], dict), "off_targets must be a dict"
    assert "near" in pr["off_targets"] and "far" in pr["off_targets"], "off_targets must have near/far"
    assert isinstance(pr["off_targets"]["near"], list), "off_targets['near'] must be a list"
    assert isinstance(pr["off_targets"]["far"], list), "off_targets['far'] must be a list"


def main():
    sequence = SEQ.replace("\n", "").replace(" ", "").upper()

    # ---------------- A ----------------
    ok, primers, info = scan_primers(sequence, size_range=(18, 25))
    assert ok, f"A failed: {info}"
    assert "forward" in primers and "reverse" in primers, "A output missing forward/reverse"

    # ---------------- B ----------------
    fwd = filter_primers(primers["forward"], min_tm=50.0, max_tm=65.0, max_hairpin_dG=-3.0, max_dimer_dG=-9.0)
    rev = filter_primers(primers["reverse"], min_tm=50.0, max_tm=65.0, max_hairpin_dG=-3.0, max_dimer_dG=-9.0)
    assert isinstance(fwd, list) and isinstance(rev, list), "B output must be lists"
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
    assert isinstance(pairs, list), "C output must be a list"
    assert len(pairs) > 0, "C produced no primer pairs"

    # 关键：只取前 2 个 pair 做联网 D（减少请求量，确保测试更容易结束）
    pairs_small = pairs[:2]

    # ---------------- D (联网) ----------------
    # 注意事项：
    # - max_workers 用 1：减少线程等待导致的“看似卡住”
    # - sleep_between_requests > 0：降低触发 NCBI 限流概率
    blast_results = process_c_module_output_d(
        c_module_output=pairs_small,
        organism="Arabidopsis thaliana",
        db="refseq_rna",
        max_workers=1,
        target_length=tl,
        near_bp=30,
        target_subject_id=None,
        min_amplicon_len=50,
        max_amplicon_len=3000,
        require_strict_orientation=True,
        sleep_between_requests=0.4,
    )

    assert isinstance(blast_results, list), "D output must be a list"
    assert len(blast_results) == len(pairs_small), "D results count must match input pairs"

    for pr in blast_results:
        assert isinstance(pr, dict), "Each PairResult must be a dict"
        assert_pairresult_contract(pr)

    print("ABCD online minimal test: PASS")
    print(f"A primers: fwd={len(primers['forward'])}, rev={len(primers['reverse'])}")
    print(f"C pairs total={len(pairs)}, tested_in_D={len(pairs_small)}")
    print("D structure contract satisfied.")


if __name__ == "__main__":
    main()
