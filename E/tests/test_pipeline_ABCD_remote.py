import unittest
import time

from pipeline_ABCD import pipeline_ABCD


class TestPipelineABCD_RemoteSoybean(unittest.TestCase):
    """
    远程（联网）测试：大豆 (Glycine max)
    - 真打 NCBI qblast
    - 会受网络/NCBI 限流影响，跑得慢是正常的
    """

    def test_pipeline_abcd_remote_soybean(self):
        sequence = (
            "ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTGTTACACAAGGTGAATCTTGGAGTGATTTTCTTCTTGACCTACAAGGC"
        ).replace("\n", "").replace(" ", "").upper()

        organism = "Glycine max"
        max_workers = 2
        sleep_between_requests = 0.2

        t0 = time.time()
        result = pipeline_ABCD(
            sequence=sequence,
            organism=organism,
            db="refseq_rna",

            dummy_blast=False,
            max_workers=max_workers,
            sleep_between_requests=sleep_between_requests,

            require_strict_orientation=True,
            min_amplicon_len=50,
            max_amplicon_len=3000,
            near_bp=30,
            
            # 保持你当前 C 阈值；如果 pairs 太少/太多你再调
            max_tm_diff=2.0,
            max_heterodimer_dG=0.0,
        )
        elapsed = time.time() - t0

        self.assertIsInstance(result, dict)
        self.assertIn("pairs", result)
        self.assertIn("blast_results", result)

        pairs = result["pairs"]
        blast_results = result["blast_results"]
        self.assertIsInstance(pairs, list)
        self.assertIsInstance(blast_results, list)
        self.assertEqual(len(blast_results), len(pairs))

        print("\n========== REMOTE RUN SUMMARY ==========")
        print(f"Organism: {organism}")
        print(f"Sequence length (target_length): {len(sequence)}")
        print(f"Pairs generated: {len(pairs)}")
        print(f"Elapsed: {elapsed:.2f}s")
        print("=======================================\n")

        if len(pairs) == 0:
            print("No primer pairs generated. Consider relaxing B/C thresholds.")
            return

        success_count = 0

        for br in blast_results:
            self.assertIn("status", br)
            self.assertIn(br["status"], ["success", "failed"])
            self.assertIn("hits", br)
            self.assertIn("amplicons", br)
            self.assertIn("target_amplicons", br)
            self.assertIn("off_targets", br)
            self.assertIn("near", br["off_targets"])
            self.assertIn("far", br["off_targets"])

            pair_idx = br.get("pair_index")
            status = br["status"]
            prod = br.get("product_size", None)

            hits_fwd = len(br["hits"]["forward"])
            hits_rev = len(br["hits"]["reverse"])

            amp_n = len(br.get("amplicons", []))
            tgt_n = len(br.get("target_amplicons", []))
            near_n = len(br.get("off_targets", {}).get("near", []))
            far_n = len(br.get("off_targets", {}).get("far", []))

            fwd = br.get("forward", "")
            rev = br.get("reverse", "")

            print(f"[pair {pair_idx}] status={status} product_size={prod}")
            print(f"  F: {fwd}")
            print(f"  R: {rev}")
            print(f"  hits_fwd={hits_fwd} hits_rev={hits_rev}")
            print(f"  amplicons={amp_n}  target_amplicons={tgt_n}  off_near={near_n}  off_far={far_n}")

            if status == "failed":
                print(f"  ERROR: {br.get('error')}\n")
                continue

            success_count += 1

            # 如果有 target_amplicons，打印第一条
            if tgt_n > 0:
                top = br["target_amplicons"][0]
                print(f"  target_subject_id={top.get('subject_id')}, gene_id={top.get('gene_id')}, length={top.get('length')}")
            elif amp_n > 0:
                top = br["amplicons"][0]
                print(f"  top_amplicon_subject_id={top.get('subject_id')}, gene_id={top.get('gene_id')}, length={top.get('length')}")

            # 分桶逻辑校验（如果 product_size 是 int）
            if isinstance(prod, int):
                for a in br["off_targets"]["near"]:
                    self.assertLessEqual(abs(a["length"] - prod), 30)
                for a in br["off_targets"]["far"]:
                    self.assertGreater(abs(a["length"] - prod), 30)

            print("")

        # 至少有一些 success（否则就是网络/限流/NCBI 解析全挂）
        self.assertGreater(success_count, 0, "All remote BLAST tasks failed; check network/NCBI throttling.")


if __name__ == "__main__":
    unittest.main(verbosity=2)
