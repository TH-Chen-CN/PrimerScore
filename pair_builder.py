from typing import List, Dict
import primer3

def build_pairs(
    forward_primers: List[Dict],
    reverse_primers: List[Dict],
    target_length: int,
    max_tm_diff: float = 3.0,  # 放宽Tm差异容忍度
    max_heterodimer_dG: float = 0.0  # 放宽异源二聚体ΔG的限制
) -> List[Dict]:
    """
    从所有引物中生成引物对，警告作为引物注释存在，并不会影响配对过程。
    
    参数：
    forward_primers (List[Dict]): 正向引物列表
    reverse_primers (List[Dict]): 反向引物列表
    target_length (int): 目标扩增片段长度
    max_tm_diff (float): 正向和反向引物Tm差异的最大容忍值
    max_heterodimer_dG (float): 异源二聚体的ΔG最大容忍值
    
    返回：
    List[Dict]: 符合条件的引物对列表
    """
    primer_pairs = []

    # 遍历所有正向引物和反向引物组合
    for fwd in forward_primers:
        for rev in reverse_primers:
            
            # 获取Tm值
            tm_fwd = fwd.get('tm_primer3', fwd.get('tm', 0))  # 设置默认值0
            tm_rev = rev.get('tm_primer3', rev.get('tm', 0))  # 设置默认值0
            tm_diff = abs(tm_fwd - tm_rev)
            
            # 打印调试信息：Tm差异
            print(f"Checking pair: Forward Tm = {tm_fwd}, Reverse Tm = {tm_rev}, Tm Diff = {tm_diff}")
            
            # 过滤：Tm差异过大
            if tm_diff > max_tm_diff:
                print(f"Skipping pair due to Tm difference: {fwd['sequence']} + {rev['sequence']} (Tm diff = {tm_diff})")
                continue
            
            # 计算异源二聚体ΔG
            hetero = primer3.calcHeterodimer(fwd['sequence'], rev['sequence'])
            hetero_dG = hetero.dg / 1000  # 除以1000转换单位为kcal/mol
            
            # 打印调试信息：异源二聚体ΔG
            print(f"Heterodimer ΔG = {hetero_dG} kcal/mol")
            
            # 过滤：异源二聚体ΔG过大
            if hetero_dG > max_heterodimer_dG:
                print(f"Skipping pair due to heterodimer ΔG: {fwd['sequence']} + {rev['sequence']} (ΔG = {hetero_dG})")
                continue
            
            # 计算扩增片段长度
            product_size = target_length

            # 打印调试信息：成功配对
            print(f"Pair accepted: {fwd['sequence']} + {rev['sequence']}")

            # 添加符合条件的引物对
            primer_pairs.append({
                "forward": fwd,
                "reverse": rev,
                "product_size": product_size,
                "tm_diff": tm_diff,
                "heterodimer_dG": hetero_dG,
                "forward_warning": fwd.get('warning', ''),  # 获取正向引物的警告
                "reverse_warning": rev.get('warning', '')   # 获取反向引物的警告
            })

    print(f"Total pairs generated: {len(primer_pairs)}")
    return primer_pairs
