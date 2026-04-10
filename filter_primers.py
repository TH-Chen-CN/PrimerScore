from typing import List, Dict
import primer3  # 导入 primer3-py 库，用于Tm计算和其他热力学计算

def filter_primers(primers: List[Dict], min_tm: float = 50.0, max_tm: float = 65.0, max_hairpin_dG: float = -3.0, max_dimer_dG: float = -9.0, min_gc: float = 0.4, max_gc: float = 0.6, min_len: int = 18, max_len: int = 25) -> List[Dict]:
    """
    对引物进行筛选，并在每个引物中添加`warning`字段，警告引物会有警告信息，合法引物没有警告。
    """
    for primer in primers:
        sequence = primer.get("sequence", "")
        primer["warning"] = ""  # 默认没有警告
        
        # 检查长度是否符合
        if len(sequence) < min_len or len(sequence) > max_len:
            primer["warning"] = "Length out of range"
            continue
        
        if 'AAAA' in sequence or 'TTTT' in sequence:
            primer["warning"] = "Too many consecutive bases"
            continue

        gc_content = primer.get("gc", 0)
        tm = primer3.calcTm(sequence)  # 使用primer3来计算Tm值
        hairpin_result = primer3.calcHairpin(sequence)  # 计算发卡结构的ΔG
        dimer_result = primer3.calcHomodimer(sequence)  # 计算自二聚体的ΔG

        # 提取ΔG值
        hairpin_dg = hairpin_result.dg  # 获取发卡结构的ΔG值
        dimer_dg = dimer_result.dg  # 获取二聚体的ΔG值

        # 检查Tm范围
        if tm < min_tm or tm > max_tm:
            primer["warning"] = "Tm out of optimal range"
            continue

        # 检查GC含量
        if gc_content < min_gc or gc_content > max_gc:
            primer["warning"] = "GC content out of optimal range"
            continue

        # 检查发卡结构
        if hairpin_dg < max_hairpin_dG:
            primer["warning"] = "Hairpin structure too strong"
            continue

        # 检查二聚体
        if dimer_dg < max_dimer_dG:
            primer["warning"] = "Dimer structure too strong"
            continue

    return primers
