import unittest

from pipeline_ABCD import pipeline_ABCD


class TestPipelineABCD_NewD(unittest.TestCase):
    def test_pipeline_abcd_dummy_blast_schema(self):
        # 输入 DNA 序列（与你之前测试一致）
        sequence = (
            "ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTG"
        ).replace("\n", "").replace(" ", "")

        # dummy_blast=True：不联网，稳定测试 schema
        result = pipeline_ABCD(
            sequence=sequence,
            organism="Arabidopsis thaliana",
            dummy_blast=True,
            max_workers=4,
            require_strict_orientation=True,
            min_amplicon_len=50,
            max_amplicon_len=3000,
            near_bp=30,
        )

        # ---- Pipeline 顶层结构 ----
        self.assertIsInstance(result, dict)
        self.assertIn("scan", result)
        self.assertIn("filter", result)
        self.assertIn("pairs", result)
        self.assertIn("blast_results", result)

        self.assertIsInstance(result["pairs"], list)
        self.assertIsInstance(result["blast_results"], list)

        # 如果 pair 为 0，D 当然也会是空 list（这不是 D 的错）
        # 但我们至少检查：blast_results 长度应与 pairs 一致（新版 pipeline 就是这样设计）
        self.assertEqual(len(result["blast_results"]), len(result["pairs"]))

        # ---- D 输出结构（新版 schema） ----
        for i, blast_result in enumerate(result["blast_results"]):
            self.assertIsInstance(blast_result, dict)

            # 基础字段
            for key in [
                "pair_index",
                "forward",
                "reverse",
                "product_size",
                "tm_diff",
                "heterodimer_dG",
                "hits",
                "amplicons",
                "target_amplicons",
                "off_targets",
                "status",
            ]:
                self.assertIn(key, blast_result)

            self.assertEqual(blast_result["pair_index"], i)
            self.assertIn(blast_result["status"], ["success", "failed"])

            # hits 结构
            self.assertIsInstance(blast_result["hits"], dict)
            self.assertIn("forward", blast_result["hits"])
            self.assertIn("reverse", blast_result["hits"])
            self.assertIsInstance(blast_result["hits"]["forward"], list)
            self.assertIsInstance(blast_result["hits"]["reverse"], list)

            # amplicons / target_amplicons / off_targets 结构
            self.assertIsInstance(blast_result["amplicons"], list)
            self.assertIsInstance(blast_result["target_amplicons"], list)
            self.assertIsInstance(blast_result["off_targets"], dict)
            self.assertIn("near", blast_result["off_targets"])
            self.assertIn("far", blast_result["off_targets"])
            self.assertIsInstance(blast_result["off_targets"]["near"], list)
            self.assertIsInstance(blast_result["off_targets"]["far"], list)

            # failed 分支：必须给 error
            if blast_result["status"] == "failed":
                self.assertIn("error", blast_result)
                self.assertIsInstance(blast_result["error"], str)
                # 失败时其它字段允许为空，但结构仍应存在（上面已检查）
                continue

            # success 分支：进一步检查 amplicon/hit 结构（dummy 模式下 hits/amplicons 不一定非空，但若存在则字段必须齐）
            # 1) 检查 hits 的字段（若有）
            for hit in blast_result["hits"]["forward"] + blast_result["hits"]["reverse"]:
                self.assertIsInstance(hit, dict)
                for hk in ["subject_id", "title", "identity", "e_value", "align_len", "sbjct_start", "sbjct_end", "strand"]:
                    self.assertIn(hk, hit)

            # 2) 检查 amplicons 的字段（若有）
            def assert_amplicon_schema(a: dict):
                self.assertIsInstance(a, dict)
                for ak in ["subject_id", "gene_id", "length", "fwd", "rev", "identity_min", "e_value_max"]:
                    self.assertIn(ak, a)
                self.assertIsInstance(a["fwd"], dict)
                self.assertIsInstance(a["rev"], dict)
                for hk in ["subject_id", "title", "identity", "e_value", "align_len", "sbjct_start", "sbjct_end", "strand"]:
                    self.assertIn(hk, a["fwd"])
                    self.assertIn(hk, a["rev"])
                self.assertIsInstance(a["length"], int)

            for a in blast_result["amplicons"]:
                assert_amplicon_schema(a)
            for a in blast_result["target_amplicons"]:
                assert_amplicon_schema(a)
            for a in blast_result["off_targets"]["near"]:
                assert_amplicon_schema(a)
            for a in blast_result["off_targets"]["far"]:
                assert_amplicon_schema(a)

            # 3) near/far 分桶逻辑（如果 product_size 存在）
            # target_length 在 pipeline 里默认=输入序列长度，所以一般会有
            tl = blast_result.get("product_size")
            if isinstance(tl, int):
                for a in blast_result["off_targets"]["near"]:
                    self.assertLessEqual(abs(a["length"] - tl), 30)
                for a in blast_result["off_targets"]["far"]:
                    self.assertGreater(abs(a["length"] - tl), 30)


if __name__ == "__main__":
    unittest.main(verbosity=2)
