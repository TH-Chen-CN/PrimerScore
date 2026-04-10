import unittest
from pipeline import run_basic_pipeline  # 导入已有的Pipeline
# 替换为并行版BLAST函数（核心修改1）
from BLAST import process_c_module_output_parallel as process_c_module_output  

class TestPipelineAndBLAST(unittest.TestCase):

    def test_pipeline_and_blast(self):
        # 输入的DNA序列
        sequence = "ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTG"

        # 调用完整Pipeline
        pipeline_result = run_basic_pipeline(sequence)
        
        # 从Pipeline的输出获取C模块结果
        c_module_output = pipeline_result['pairs']

        # 调用并行版BLAST模块（核心修改2：新增max_workers参数，控制并发）
        organism = "Arabidopsis thaliana"  # 用户输入的物种
        blast_results = process_c_module_output(
            c_module_output=c_module_output,
            organism=organism,
            max_workers=5  # 测试环境建议降低并发数，避免NCBI限流
        )

        # 验证BLAST结果基础格式
        self.assertIsInstance(blast_results, list)
        
        # 验证每对引物的数据
        for result in blast_results:
            # 基础字段校验（兼容并行版返回格式）
            self.assertIn("forward", result)
            self.assertIn("reverse", result)
            self.assertIn("target_amplicons", result)  # 核心修改3：字段名从target_amplicon→target_amplicons
            self.assertIn("status", result)  # 新增：校验任务状态（success/failed）

            # 状态分支校验：成功场景
            if result["status"] == "success":
                # 检查扩增片段字段类型
                self.assertIsInstance(result["target_amplicons"], list)

                # 检查目标扩增片段的核心字段（若有结果）
                if result["target_amplicons"]:
                    amplicon = result["target_amplicons"][0]
                    self.assertIn("gene_id", amplicon)
                    self.assertIn("align_len", amplicon)  # 核心修改4：字段名从length→align_len（并行版命名）
                    self.assertIn("strand", amplicon)
                    self.assertIn("identity", amplicon)  # 新增：校验一致性百分比
                    self.assertIn("e_value", amplicon)   # 新增：校验e值（显著性）
            
            # 状态分支校验：失败场景（容错性校验）
            elif result["status"] == "failed":
                self.assertIn("error", result)  # 失败时必须包含错误信息
                self.assertIsInstance(result["error"], str)  # 错误信息为字符串类型

            # 非法状态校验（防御性编程）
            else:
                self.fail(f"无效的BLAST结果状态：{result['status']}")

if __name__ == '__main__':
    # 运行测试（新增verbosity=2，输出详细测试日志）
    unittest.main(verbosity=2)