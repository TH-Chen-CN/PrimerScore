import primer3  # 导入 primer3-py 库
from typing import Tuple, List, Dict

def revcomp(seq: str) -> str:
    """计算互补序列（不反转，只是计算互补）"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in seq)  # 直接取互补，不反转

def scan_primers(sequence: str, size_range: Tuple[int, int] = (18, 25)) -> Tuple[bool, Dict, Dict]:
    """扫描并生成引物"""
    
    # 输入验证
    if not isinstance(sequence, str) or len(sequence) < 100 or not all(base in 'ACGT' for base in sequence):
        return False, {}, {"error": "Invalid sequence. Must be at least 100 bases long and contain only ACGT."}

    min_len, max_len = size_range
    forward_primers = []
    reverse_primers = []

    # 生成正向引物（从 5' 端开始）
    for L in range(min_len, max_len + 1):
        primer = sequence[:L]  # 直接从 5' 端起始生成引物
        gc_content = (primer.count('G') + primer.count('C')) / L
        tm = primer3.calcTm(primer)  # 使用 primer3 计算 Tm
        forward_primers.append({
            "sequence": primer,
            "start": 0,
            "end": L,
            "length": L,
            "gc": gc_content,
            "tm": tm
        })
    
    # 生成反向引物（从 3' 端开始，直接互补并反向）
    for L in range(min_len, max_len + 1):
        # 从 3'端起始，截取 primer_length 长度的序列
        reverse_primer_sequence = sequence[-L:]  # 取目标序列的末尾部分
        reverse_primer = revcomp(reverse_primer_sequence[::-1])  # 先反向再互补
        
        gc_content = (reverse_primer.count('G') + reverse_primer.count('C')) / L
        tm = primer3.calcTm(reverse_primer)  # 使用 primer3 计算 Tm
        reverse_primers.append({
            "sequence": reverse_primer,
            "start": len(sequence) - L,
            "end": len(sequence),
            "length": L,
            "gc": gc_content,
            "tm": tm
        })

    return True, {"forward": forward_primers, "reverse": reverse_primers}, {"info": "Primer generation successful"}
