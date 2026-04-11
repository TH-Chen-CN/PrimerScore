# PrimerScore

**PCR primer design with BLAST-based specificity checking, expression-aware evaluation, and combined scoring.**

PCR 引物设计工具，集成 BLAST 特异性分析、表达量感知评估与综合打分排序。

> **Current species support | 当前物种支持：Soybean (*Glycine max*) only**

---

## Overview | 简介

PrimerScore is a desktop application that automates the full PCR primer design workflow:

1. Scans a target DNA sequence for candidate primers
2. Filters by thermodynamic criteria (Tm, GC%, hairpin, dimer)
3. Assembles valid primer pairs
4. Checks specificity via **remote NCBI BLAST** (internet required)
5. Evaluates each pair against gene expression profiles (SoyOmics + CoNekT)
6. Produces a ranked table of primer pairs with warnings and export options

PrimerScore 是一个桌面工具，自动完成从候选引物扫描到综合打分排序的完整 PCR 引物设计流程。BLAST 分析通过联网调用 **NCBI 远程 BLAST 服务**完成，无需本地数据库。

---

## Important: No Local Database Required | 重要：无需本地数据库

PrimerScore uses **NCBI's public BLAST web service** for specificity analysis. There is **no local sequence database** bundled with this tool. Every run that uses real BLAST mode will make network requests to NCBI servers.

本工具使用 **NCBI 公开 BLAST 网络服务**进行特异性分析，**不包含本地序列数据库**。真实 BLAST 模式下每次运行均需联网调用 NCBI 接口。

**For researchers who want faster or offline analysis | 希望进行本地部署的研究者：**

The BLAST module (`blast_client_d.py`) is designed with a clean interface that can be adapted to query a local BLAST+ installation or institution-specific database. Researchers are welcome to fork this project and integrate their own local BLAST databases — contributions and adaptations are encouraged.

BLAST 模块（`blast_client_d.py`）接口清晰，可改造为调用本地 BLAST+ 安装或机构自有数据库。欢迎有意向的专家学者 fork 本项目，将其移植到本地数据库环境中使用，期待您的改进与贡献。

---

## Features | 功能

- **Fixed-end primer scanning** — scans across the full input sequence for forward/reverse candidates (18–25 bp)  
  固定端引物扫描，遍历全序列生成候选正向/反向引物

- **Thermodynamic filtering** — Tm range (50–65°C), GC content, hairpin and homodimer checks  
  热力学筛选：Tm 范围、GC 含量、发夹结构与同源二聚体检测

- **Primer pair assembly** — Tm difference and heterodimer screening  
  引物对构建，包含 Tm 差值与异源二聚体筛查

- **BLAST specificity analysis** — remote NCBI BLAST against RefSeq RNA, off-target detection  
  远程 NCBI BLAST 特异性分析，检测脱靶情况

- **Expression-aware evaluation** — integrates SoyOmics and CoNekT expression data across tissues  
  表达量感知评估，整合 SoyOmics 与 CoNekT 多组织表达数据

- **Combined ranking** — scores primer pairs by thermodynamics, expression, and specificity  
  综合打分排序，综合热力学、表达量与特异性指标

- **Desktop GUI** with result table, double-click detail view, copyable sequences, JSON/CSV export  
  桌面图形界面，支持结果表格展示、双击查看详情、序列复制、JSON/CSV 导出

---

## Installation | 安装

### Requirements | 依赖

- Python 3.9+
- Internet access (for BLAST and expression queries) | 联网（用于 BLAST 和表达量查询）

### Steps | 安装步骤

```bash
# 1. Clone the repository | 克隆仓库
git clone https://github.com/TH-Chen-CN/PrimerScore.git
cd PrimerScore

# 2. (Recommended) Create a virtual environment | 建议使用虚拟环境
python -m venv venv
source venv/bin/activate        # Linux/macOS
venv\Scripts\activate           # Windows

# 3. Install dependencies | 安装依赖
pip install -r requirements.txt
```

**Linux users**: tkinter may need to be installed separately:
```bash
sudo apt-get install python3-tk   # Debian/Ubuntu
sudo dnf install python3-tkinter  # Fedora/RHEL
```

---

## Usage | 使用方法

### Launch the GUI | 启动图形界面

```bash
python gui.py
```

### Step-by-step | 操作步骤

1. **Paste your DNA sequence** into the "Target Sequence" box.  
   Raw sequence or FASTA format both accepted. Non-ACGT characters are stripped automatically.  
   将目标 DNA 序列粘贴到"Target Sequence"框中，支持原始序列或 FASTA 格式，非 ACGT 字符自动去除。

2. **Set the organism** (default: *Glycine max*).  
   设置物种（默认：大豆 *Glycine max*）。

3. **Set Top N** — how many top-ranked primer pairs to display (default: 10).  
   设置 Top N，即展示的最优引物对数量（默认：10）。

4. **Choose BLAST mode**:
   - **Dummy mode (checked)** — skips BLAST, fast, for testing the pipeline only.  
     Dummy 模式（勾选）：跳过 BLAST，速度快，仅用于流程测试。
   - **Real mode (unchecked)** — runs full remote BLAST. Requires internet. Takes 1–5 minutes.  
     真实模式（不勾选）：运行完整远程 BLAST，需联网，耗时约 1–5 分钟。

5. **Click "Run PrimerScore"** to start. Progress is shown in the Log Output panel.  
   点击"Run PrimerScore"开始运行，进度显示在 Log Output 面板。

6. **Review results** in the Combined Ranking Results table. Double-click any row for full details.  
   在结果表格中查看排名，双击任意行查看详情。

7. **Export** results using "Export JSON" or "Export CSV".  
   使用"Export JSON"或"Export CSV"导出结果。

> **Tip | 提示**: Always test with Dummy mode first to verify your sequence is parsed correctly, then uncheck it for a real BLAST run.  
> 建议先用 Dummy 模式验证流程，确认无误后再取消勾选进行真实 BLAST 分析。

---

## Pipeline Architecture | 流程架构

PrimerScore runs a 6-stage pipeline (ABCDES):

| Stage | Module | Description |
|-------|--------|-------------|
| **A** | `scan_primers.py` | Scan target sequence for candidate primers (18–25 bp) |
| **B** | `filter_primers.py` | Thermodynamic filtering (Tm, GC%, hairpin, homodimer) |
| **C** | `pair_builder.py` | Pair assembly with Tm-delta and heterodimer checks |
| **D** | `blast_client_d.py` | Remote NCBI BLAST for off-target and specificity analysis |
| **E** | `E/evaluate.py` | Expression-aware evaluation (SoyOmics + CoNekT integration) |
| **S** | `scoring.py` | Combined ranking: thermo score + expression score |

The expression module (E) resolves BLAST hit IDs (RefSeq XM_/NM_) to Glyma gene models via NCBI E-utilities, then queries expression databases. Results are cached locally in SQLite to reduce repeated network calls.

流程分为 6 个阶段（ABCDES），表达量模块（E）会将 BLAST 命中 ID（RefSeq XM_/NM_）解析为 Glyma 基因模型，再查询表达数据库，结果通过本地 SQLite 缓存以减少重复网络请求。

---

## Project Structure | 项目结构

```text
PrimerScore/
├── gui.py                      # GUI entry point | 图形界面入口
├── pipeline_ABCDES.py          # Pipeline orchestrator | 流程总调度
├── scan_primers.py             # Module A: primer scanning | 模块A：引物扫描
├── filter_primers.py           # Module B: thermodynamic filter | 模块B：热力学筛选
├── pair_builder.py             # Module C: pair assembly | 模块C：引物对构建
├── blast_client_d.py           # Module D: BLAST specificity | 模块D：BLAST特异性
├── scoring.py                  # Module S: ranking & scoring | 模块S：打分排序
├── PrimerScore.spec            # PyInstaller build config | 打包配置
├── requirements.txt            # Python dependencies | Python依赖
└── E/                          # Expression evaluation module | 表达量评估模块
    ├── evaluate.py             # Module E: expression evaluation | 模块E：表达量评估
    ├── expression/
    │   ├── client.py           # Unified expression client | 统一表达量客户端
    │   ├── soyomics_client.py  # SoyOmics data fetcher | SoyOmics 数据获取
    │   ├── conekt_client.py    # CoNekT data fetcher | CoNekT 数据获取
    │   ├── integrator.py       # Multi-source data merging | 多源数据融合
    │   ├── canonical.py        # Tissue condition normalization | 组织条件规范化
    │   ├── cache.py            # SQLite caching layer | SQLite 缓存层
    │   ├── atlas.py            # Expression atlas utilities | 表达图谱工具
    │   └── fusion.py           # Expression data fusion | 表达数据融合
    └── resolver/
        ├── resolver.py         # RefSeq → Glyma ID resolver | ID解析器
        └── cache.py            # Resolver SQLite cache | 解析缓存
```

---

## Building an Executable | 打包为可执行文件

PrimerScore can be packaged into a standalone executable using PyInstaller:

```bash
pip install pyinstaller
pyinstaller PrimerScore.spec
```

The output will be in the `dist/` folder. Note that executables are **platform-specific** — a Windows build must be run on Windows, a Linux build on Linux.

打包产物在 `dist/` 目录下。注意可执行文件**平台相关**，Windows 打包需在 Windows 上进行，Linux 打包需在 Linux 上进行。

---

## Notes | 注意事项

- Real BLAST mode requires a stable internet connection and may take several minutes per run.  
  真实 BLAST 模式需要稳定的网络连接，每次运行可能需要数分钟。

- Expression data is cached locally in `expression_cache.sqlite` and `resolver_cache.sqlite`. Deleting these files clears the cache and forces fresh network queries.  
  表达量数据会缓存在本地 SQLite 文件中，删除缓存文件可强制重新联网查询。

- Currently only soybean (*Glycine max*) is supported for expression-aware scoring. Other organisms will complete BLAST analysis but skip expression evaluation.  
  当前仅支持大豆（*Glycine max*）的表达量评估；其他物种可完成 BLAST 分析，但跳过表达量评分。

- If PyInstaller packaging fails due to missing imports, check `PrimerScore.spec` and add the missing module to `hiddenimports`.  
  若打包失败提示缺少模块，请在 `PrimerScore.spec` 的 `hiddenimports` 中添加对应模块名。

---

## Contributing & Local Deployment | 贡献与本地部署

This project was built to use public web APIs (NCBI BLAST, SoyOmics, CoNekT). Researchers and bioinformaticians who wish to:

- Adapt this tool to **other plant species**
- Integrate a **local BLAST+ database** for faster, offline analysis
- Add expression data from other sources

...are warmly welcomed to fork and extend this project. The modular pipeline design (ABCDES) makes it straightforward to swap out individual components.

本项目基于公开 Web API 构建。欢迎有意向的专家学者：

- 将本工具扩展至**其他植物物种**
- 集成**本地 BLAST+ 数据库**以实现更快的离线分析
- 接入其他来源的表达量数据

欢迎 fork 本项目并进行改进，模块化的 ABCDES 流程设计使各组件的替换与扩展十分方便。

---

## Contact | 联系方式

**Author | 作者**: Chan  
**Email**: cth0619hhh@gmail.com

Issues and pull requests are welcome. | 欢迎提交 Issue 和 Pull Request。
