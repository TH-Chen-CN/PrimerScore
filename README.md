# PrimerScore

PrimerScore is a desktop tool for PCR primer design with BLAST specificity checking, expression-aware evaluation, and ranking.

PrimerScore 是一个用于 PCR 引物设计的桌面工具，集成了 BLAST 特异性分析、表达量感知评估以及综合评分排序功能。

---

## 🚀 Features | 功能

- Fixed-end primer scanning  
  固定端引物扫描

- Thermodynamic filtering  
  引物热力学筛选

- Primer pair assembly  
  引物对构建

- BLAST-based specificity analysis  
  基于 BLAST 的特异性分析

- Expression-aware evaluation  
  表达量感知评估

- Combined ranking of primer pairs  
  引物对综合评分排序

- Desktop GUI with:
  桌面图形界面功能：

  - Result table  
    结果表格展示

  - Double-click detail view  
    双击查看详细信息

  - Copyable primer sequences  
    引物序列可复制

  - JSON / CSV export  
    支持 JSON / CSV 导出

---

## 📁 Project Structure | 项目结构

```text
gui.py                  GUI入口
pipeline_ABCDES.py      主流程
scan_primers.py         模块A
filter_primers.py       模块B
pair_builder.py         模块C
blast_client_d.py       模块D
E/                      表达量评估模块
scoring.py              打分模块
tests/                  测试文件
