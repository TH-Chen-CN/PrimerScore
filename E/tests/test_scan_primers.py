from scan_primers import scan_primers# 假设 scan_primers 已经从模块A导入

def test_primer_generation():
    # 输入给定的序列
    sequence = "ATGGAAGGATCATCAGGTGTAAGGAAAGGCGCATGGAGTCAAATTGAAGATAACCTTCTCAGAGATTGCGTGAACCTTCATGGGGAAGGAAAATGGCACCTTGTTCCTAAAAGAGCAGGGTTGAACAGATGCCGCAAGAGTTGTAGATTGAGATGGTTGAACTATCTTAAACCAAATATCAAGCGGGGAGATTTTAGTGAAGATGAGGTTGATTTGATGATCAGATTGCACAAGCTTTTGGGAAACAGATGGTCCCTAATTGCAGGGAGACTTCCAGGAAGAACTTCAAACGATGTGAAGAATTACTGGAACACCTACATGCGCCGTAAAGTACACTCTCACAAGAAAGACAACAACATAGAGAAGCAAGCTGATGAGGCCAAACCAATAGTGAAACATCACGAAGTTATAAAACCTGTTCCTCGAACTCTATCAAAAACATCTCCATGGTTGCAAGGGAAATTTGTTAATAGTTCCAAAGTTGGTGTTAGTGAAGAAGGTGCAACTTCAATATCAGGGTCTGCTGGGAATTGGTGGGAAACTTTGTTAGATGACAAAGAAGACAATGCAGTTAATAACAACAACACGTGCTTCTTTGGTGGGGCAGATGGAGAGTTTAACCTTTGGAGTGAAGAGCTTACTTCAATTGATTGTGATTTTGTTACACAAGGTGAATCTTGGAGTGATTTTCTTCTTGACCTACAAGGC"
    
    # 设置引物长度范围
    size_range = (18, 25)  # 设置正向和反向引物的长度范围
    
    # 调用 scan_primers 函数生成引物
    success, result, status = scan_primers(sequence, size_range)
    
    # 检查是否生成了引物
    assert success, "引物生成失败"
    print(f"生成的正向引物：{len(result['forward'])} 条")
    print(f"生成的反向引物：{len(result['reverse'])} 条")
    # 打印全部生成的引物
    print("\n全部正向引物：")
    for primer in result["forward"]:
        print(f"序列: {primer['sequence']}, 起始位置: {primer['start']}, 结束位置: {primer['end']}, GC含量: {primer['gc']}, Tm: {primer['tm']}")

    print("\n全部反向引物：")
    for primer in result["reverse"]:
        print(f"序列: {primer['sequence']}, 起始位置: {primer['start']}, 结束位置: {primer['end']}, GC含量: {primer['gc']}, Tm: {primer['tm']}")
    # 输出状态信息
    print("\n状态信息：", status)

# 运行测试
test_primer_generation()
