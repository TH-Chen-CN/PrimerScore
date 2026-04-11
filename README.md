# PrimerScore

**PCR primer design with BLAST-based specificity checking, expression-aware evaluation, and combined scoring.**

[中文版本](#中文版本)

> **Current species support: Soybean (*Glycine max*) only**

---

## Overview

PrimerScore is a desktop application that automates the full PCR primer design workflow:

1. Scans a target DNA sequence for candidate primers
2. Filters by thermodynamic criteria (Tm, GC%, hairpin, dimer)
3. Assembles valid primer pairs
4. Checks specificity via **remote NCBI BLAST** (internet required)
5. Evaluates each pair against gene expression profiles (SoyOmics + CoNekT)
6. Produces a ranked table of primer pairs with warnings and export options

---

## Important: No Local Database Required

PrimerScore uses **NCBI's public BLAST web service** for specificity analysis. There is **no local sequence database** bundled with this tool. Every run that uses real BLAST mode will make network requests to NCBI servers.

**For researchers who want faster or offline analysis:**

The BLAST module (`blast_client_d.py`) is designed with a clean interface that can be adapted to query a local BLAST+ installation or institution-specific database. Researchers are welcome to fork this project and integrate their own local BLAST databases — contributions and adaptations are encouraged.

---

## Features

- **Fixed-end primer scanning** — scans across the full input sequence for forward/reverse candidates (18–25 bp)
- **Thermodynamic filtering** — Tm range (50–65°C), GC content, hairpin and homodimer checks
- **Primer pair assembly** — Tm difference and heterodimer screening
- **BLAST specificity analysis** — remote NCBI BLAST against RefSeq RNA, off-target detection
- **Expression-aware evaluation** — integrates SoyOmics and CoNekT expression data across tissues
- **Combined ranking** — scores primer pairs by thermodynamics, expression, and specificity
- **Desktop GUI** with result table, double-click detail view, copyable sequences, JSON/CSV export

---

## Installation

### Requirements

- Python 3.9+
- Internet access (for BLAST and expression queries)

### Steps

```bash
# 1. Clone the repository
git clone https://github.com/TH-Chen-CN/PrimerScore.git
cd PrimerScore

# 2. (Recommended) Create a virtual environment
python -m venv venv
source venv/bin/activate        # Linux/macOS
venv\Scripts\activate           # Windows

# 3. Install dependencies
pip install -r requirements.txt
```

**Linux users**: tkinter may need to be installed separately:
```bash
sudo apt-get install python3-tk   # Debian/Ubuntu
sudo dnf install python3-tkinter  # Fedora/RHEL
```

---

## Usage

### Launch the GUI

```bash
python gui.py
```

### Step-by-step

1. **Paste your DNA sequence** into the "Target Sequence" box.
   Raw sequence or FASTA format both accepted. Non-ACGT characters are stripped automatically.

2. **Set the organism** (default: *Glycine max*).

3. **Set Top N** — how many top-ranked primer pairs to display (default: 10).

4. **Choose BLAST mode**:
   - **Dummy mode (checked)** — skips BLAST, fast, for testing the pipeline only.
   - **Real mode (unchecked)** — runs full remote BLAST. Requires internet. Takes 1–5 minutes.

5. **Click "Run PrimerScore"** to start. Progress is shown in the Log Output panel.

6. **Review results** in the Combined Ranking Results table. Double-click any row for full details.

7. **Export** results using "Export JSON" or "Export CSV".

> **Tip**: Always test with Dummy mode first to verify your sequence is parsed correctly, then uncheck it for a real BLAST run.

---

## Pipeline Architecture

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

---

## Project Structure

```text
PrimerScore/
├── gui.py                      # GUI entry point
├── pipeline_ABCDES.py          # Pipeline orchestrator
├── scan_primers.py             # Module A: primer scanning
├── filter_primers.py           # Module B: thermodynamic filter
├── pair_builder.py             # Module C: pair assembly
├── blast_client_d.py           # Module D: BLAST specificity
├── scoring.py                  # Module S: ranking & scoring
├── PrimerScore.spec            # PyInstaller build config
├── requirements.txt            # Python dependencies
└── E/                          # Expression evaluation module
    ├── evaluate.py             # Module E: expression evaluation
    ├── expression/
    │   ├── client.py           # Unified expression client
    │   ├── soyomics_client.py  # SoyOmics data fetcher
    │   ├── conekt_client.py    # CoNekT data fetcher
    │   ├── integrator.py       # Multi-source data merging
    │   ├── canonical.py        # Tissue condition normalization
    │   ├── cache.py            # SQLite caching layer
    │   ├── atlas.py            # Expression atlas utilities
    │   └── fusion.py           # Expression data fusion
    └── resolver/
        ├── resolver.py         # RefSeq → Glyma ID resolver
        └── cache.py            # Resolver SQLite cache
```

---

## Building an Executable

PrimerScore can be packaged into a standalone executable using PyInstaller:

```bash
pip install pyinstaller
pyinstaller PrimerScore.spec
```

The output will be in the `dist/` folder. Note that executables are **platform-specific** — a Windows build must be run on Windows, a Linux build on Linux.

---

## Notes

- Real BLAST mode requires a stable internet connection and may take several minutes per run.
- Expression data is cached locally in `expression_cache.sqlite` and `resolver_cache.sqlite`. Deleting these files clears the cache and forces fresh network queries.
- Currently only soybean (*Glycine max*) is supported for expression-aware scoring. Other organisms will complete BLAST analysis but skip expression evaluation.
- If PyInstaller packaging fails due to missing imports, check `PrimerScore.spec` and add the missing module to `hiddenimports`.

---

## Contributing & Local Deployment

This project was built to use public web APIs (NCBI BLAST, SoyOmics, CoNekT). Researchers and bioinformaticians who wish to:

- Adapt this tool to **other plant species**
- Integrate a **local BLAST+ database** for faster, offline analysis
- Add expression data from other sources

...are warmly welcomed to fork and extend this project. The modular pipeline design (ABCDES) makes it straightforward to swap out individual components.

---

## Contact

**Author**: Chan
**Email**: cth0619hhh@gmail.com

Issues and pull requests are welcome.

---

---

# 中文版本

**PCR 引物设计工具，集成 BLAST 特异性分析、表达量感知评估与综合打分排序。**

[English Version](#primerscore)

> **当前物种支持：仅大豆（*Glycine max*）**

---

## 简介

PrimerScore 是一个桌面工具，自动完成从候选引物扫描到综合打分排序的完整 PCR 引物设计流程：

1. 扫描目标 DNA 序列，生成候选引物
2. 按热力学标准筛选（Tm、GC%、发夹结构、二聚体）
3. 构建有效引物对
4. 通过**远程 NCBI BLAST** 检查特异性（需联网）
5. 结合基因表达谱评估每对引物（SoyOmics + CoNekT）
6. 输出带警告信息和导出功能的引物对排名表

---

## 重要：无需本地数据库

PrimerScore 使用 **NCBI 公开 BLAST 网络服务**进行特异性分析，**不包含本地序列数据库**。真实 BLAST 模式下每次运行均需联网调用 NCBI 接口。

**希望进行本地部署的研究者：**

BLAST 模块（`blast_client_d.py`）接口清晰，可改造为调用本地 BLAST+ 安装或机构自有数据库。欢迎有意向的专家学者 fork 本项目，将其移植到本地数据库环境中使用，期待您的改进与贡献。

---

## 功能

- **固定端引物扫描** — 遍历全序列生成候选正向/反向引物（18–25 bp）
- **热力学筛选** — Tm 范围（50–65°C）、GC 含量、发夹结构与同源二聚体检测
- **引物对构建** — 包含 Tm 差值与异源二聚体筛查
- **BLAST 特异性分析** — 远程 NCBI BLAST 比对 RefSeq RNA，检测脱靶情况
- **表达量感知评估** — 整合 SoyOmics 与 CoNekT 多组织表达数据
- **综合打分排序** — 综合热力学、表达量与特异性指标
- **桌面图形界面** — 支持结果表格展示、双击查看详情、序列复制、JSON/CSV 导出

---

## 安装

### 依赖

- Python 3.9+
- 联网（用于 BLAST 和表达量查询）

### 安装步骤

```bash
# 1. 克隆仓库
git clone https://github.com/TH-Chen-CN/PrimerScore.git
cd PrimerScore

# 2. 建议使用虚拟环境
python -m venv venv
source venv/bin/activate        # Linux/macOS
venv\Scripts\activate           # Windows

# 3. 安装依赖
pip install -r requirements.txt
```

**Linux 用户**：tkinter 可能需要单独安装：
```bash
sudo apt-get install python3-tk   # Debian/Ubuntu
sudo dnf install python3-tkinter  # Fedora/RHEL
```

---

## 使用方法

### 启动图形界面

```bash
python gui.py
```

### 操作步骤

1. 将目标 DNA 序列粘贴到 **"Target Sequence"** 框中，支持原始序列或 FASTA 格式，非 ACGT 字符自动去除。

2. **设置物种**（默认：大豆 *Glycine max*）。

3. **设置 Top N**，即展示的最优引物对数量（默认：10）。

4. **选择 BLAST 模式**：
   - **Dummy 模式（勾选）**：跳过 BLAST，速度快，仅用于流程测试。
   - **真实模式（不勾选）**：运行完整远程 BLAST，需联网，耗时约 1–5 分钟。

5. 点击 **"Run PrimerScore"** 开始运行，进度显示在 Log Output 面板。

6. 在结果表格中查看排名，双击任意行查看详情。

7. 使用 **"Export JSON"** 或 **"Export CSV"** 导出结果。

> **提示**：建议先用 Dummy 模式验证流程，确认无误后再取消勾选进行真实 BLAST 分析。

---

## 流程架构

PrimerScore 运行 6 阶段流程（ABCDES）：

| 阶段 | 模块 | 说明 |
|------|------|------|
| **A** | `scan_primers.py` | 扫描目标序列，生成候选引物（18–25 bp） |
| **B** | `filter_primers.py` | 热力学筛选（Tm、GC%、发夹结构、同源二聚体） |
| **C** | `pair_builder.py` | 引物对构建，含 Tm 差值与异源二聚体筛查 |
| **D** | `blast_client_d.py` | 远程 NCBI BLAST 脱靶与特异性分析 |
| **E** | `E/evaluate.py` | 表达量感知评估（SoyOmics + CoNekT 整合） |
| **S** | `scoring.py` | 综合排序：热力学得分 + 表达量得分 |

表达量模块（E）会将 BLAST 命中 ID（RefSeq XM_/NM_）解析为 Glyma 基因模型，再查询表达数据库，结果通过本地 SQLite 缓存以减少重复网络请求。

---

## 项目结构

```text
PrimerScore/
├── gui.py                      # 图形界面入口
├── pipeline_ABCDES.py          # 流程总调度
├── scan_primers.py             # 模块A：引物扫描
├── filter_primers.py           # 模块B：热力学筛选
├── pair_builder.py             # 模块C：引物对构建
├── blast_client_d.py           # 模块D：BLAST特异性
├── scoring.py                  # 模块S：打分排序
├── PrimerScore.spec            # 打包配置
├── requirements.txt            # Python依赖
└── E/                          # 表达量评估模块
    ├── evaluate.py             # 模块E：表达量评估
    ├── expression/
    │   ├── client.py           # 统一表达量客户端
    │   ├── soyomics_client.py  # SoyOmics 数据获取
    │   ├── conekt_client.py    # CoNekT 数据获取
    │   ├── integrator.py       # 多源数据融合
    │   ├── canonical.py        # 组织条件规范化
    │   ├── cache.py            # SQLite 缓存层
    │   ├── atlas.py            # 表达图谱工具
    │   └── fusion.py           # 表达数据融合
    └── resolver/
        ├── resolver.py         # RefSeq → Glyma ID 解析器
        └── cache.py            # 解析缓存
```

---

## 打包为可执行文件

可使用 PyInstaller 将 PrimerScore 打包为独立可执行文件：

```bash
pip install pyinstaller
pyinstaller PrimerScore.spec
```

打包产物在 `dist/` 目录下。注意可执行文件**平台相关**，Windows 打包需在 Windows 上进行，Linux 打包需在 Linux 上进行。

---

## 注意事项

- 真实 BLAST 模式需要稳定的网络连接，每次运行可能需要数分钟。
- 表达量数据会缓存在本地 `expression_cache.sqlite` 与 `resolver_cache.sqlite` 中，删除缓存文件可强制重新联网查询。
- 当前仅支持大豆（*Glycine max*）的表达量评估；其他物种可完成 BLAST 分析，但跳过表达量评分。
- 若打包失败提示缺少模块，请在 `PrimerScore.spec` 的 `hiddenimports` 中添加对应模块名。

---

## 贡献与本地部署

本项目基于公开 Web API 构建（NCBI BLAST、SoyOmics、CoNekT）。欢迎有意向的专家学者：

- 将本工具扩展至**其他植物物种**
- 集成**本地 BLAST+ 数据库**以实现更快的离线分析
- 接入其他来源的表达量数据

欢迎 fork 本项目并进行改进，模块化的 ABCDES 流程设计使各组件的替换与扩展十分方便。

---

## 联系方式

**作者**：Chan
**邮箱**：cth0619hhh@gmail.com

欢迎提交 Issue 和 Pull Request。
