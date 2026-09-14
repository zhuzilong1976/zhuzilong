# 换到另一台机器：操作指南

日期：2026-09-14

---

## 0. 一句话流程

**打包 → 拷过去 → 装 R 包 → 跑稳定性复现。** 已完成的分析结果都在包里，不需要重跑。

---

## 1. 先决定搬什么

| 内容 | 体积 | 是否必搬 |
|---|---:|---|
| `outputs/repo/`（代码、文档、图表、表格） | 1 MB | **必搬** |
| `work/data/`（counts.mtx 等派生数据） | 212 MB | **必搬**（否则要从 h5ad 重新生成，约 30 秒，但需要 h5ad） |
| 主分析结果 `outputs/results/vko_ms_microglia` | 24 MB | **建议搬**（含 4 个 800×800 全基因扰动矩阵，重建需 2 小时/网络） |
| 技术对照网络 `outputs/results/technical_null` | 11 MB | **建议搬**（重建需 2 小时/网络） |
| 参数敏感性 `outputs/results/param_sensitivity` | 5.6 MB | 建议搬 |
| 跨数据集复现 `outputs/results/replication` | < 1 MB | 建议搬 |
| 深度匹配结果 `outputs/results/vko_ms_depthmatched` | 2 MB | 建议搬（如果本机这轮完成了） |
| `work/stability/`（已完成的稳定性重复） | 视进度 | **建议搬**（新机器会跳过已完成的部分） |
| CellOracle 结果 `outputs/results/celloracle` | 282 MB | 可不搬（12 分钟可重建，但需要 Python 环境） |
| 作者网络 `work/author_networks` | 108 MB | 可不搬（2 分钟可重新下载） |
| 原始 h5ad `work/geo` | 128 MB | 可不搬（可重新下载） |

**最小必搬约 250 MB；加上建议项约 660 MB；全量约 1.1 GB。**

---

## 2. 打包（在当前机器上执行）

仓库里已写好打包脚本：

```powershell
# 常规包（约 250 MB）
powershell -File outputs\repo\code\00_setup\make_transfer_bundle.ps1 -Dest D:\transfer

# 全量包（含 CellOracle 结果与作者网络，约 1.1 GB）
powershell -File outputs\repo\code\00_setup\make_transfer_bundle.ps1 -Dest D:\transfer -Full

# 压缩
Compress-Archive -Path D:\transfer\vko-project -DestinationPath D:\transfer\vko-project.zip
```

脚本会在 `D:\transfer\vko-project\` 下生成这样的结构——**这正是脚本内部的相对路径所要求的布局**（`work/` 与 `outputs/` 与 `code/` 同级）：

```
vko-project/
├── code/            ← 全部分析脚本
├── docs/            ← 分析文档（含 Methods 初稿）
├── env/             ← 环境记录
├── results/         ← 图表与补充表
├── data/            ← 数据说明
├── work/            ← 派生数据与稳定性结果
├── outputs/results/ ← 已完成的分析输出
└── README.md / LICENSE / CITATION.cff
```

> **不要只拷 `outputs/repo/`**：脚本里的路径是 `work/data/...` 与 `outputs/results/...`（这是它们实际执行时的路径，我刻意没有改写，以保证"发布的代码 = 跑出结果的代码"）。只拷 repo 会导致所有脚本找不到输入。

---

## 3. 在新机器上装环境

新机器通常能正常联网，**比本机简单得多**（本机当初是因为沙箱无法用 TLS，才走了下载二进制包离线安装的绕路）。

### 3.1 R 与依赖

要求 R ≥ 4.4（本机用的是 4.6.1）。

```r
install.packages("scTenifoldKnk")          # CRAN，v1.1+
install.packages(c("scTenifoldNet", "Matrix", "MASS", "cli", "igraph",
                   "reshape2", "enrichR", "locfdr", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("fgsea", "preprocessCore"))
```

验证（这一步会跑一次小规模端到端测试）：

```bash
Rscript code/00_setup/00_setup_environment.R
```

### 3.2 Python（仅当要重跑 CellOracle 对照时需要）

在 Linux 上 `pip install celloracle` 通常可以直接成功（`velocyto`、`gimmemotifs` 在 Linux 有源码可编译，需要 `gcc`、`gsl` 等）。本机在 Windows 上无法编译，才用了替身模块。

**如果只想补稳定性复现，完全不需要装 Python。**

---

## 4. 在新机器上跑稳定性复现

这是换机器的主要目标。仓库里已提供三种跑法：

### Linux / macOS（推荐）

```bash
cd vko-project
chmod +x code/03_controls/run_stability_linux.sh
./code/03_controls/run_stability_linux.sh 8 10 800 800
#                                    ↑  ↑   ↑   ↑
#                          并行数 ┘  │   │   └ 每重复细胞数
#                            重复数 ─┘   └ 网络基因数
```

### 有 Slurm/PBS 调度器

```bash
mkdir -p logs && sbatch code/03_controls/run_stability_slurm.sh
```

### Windows

```powershell
powershell -File code/03_controls/run_stability_local.ps1 -Parallel 4 -Replicates 10 -Genes 800 -Rscript "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
```

### 汇总（三种跑法都需要）

```bash
Rscript code/03_controls/29_stability_analysis.R
```

输出：

- `results/tables/Table_S10_stability_replicates.csv` — 逐重复的合规率、网络密度、出度 0 基因数
- `results/tables/Table_S11_stability_summary.csv` — 四个检查点的判定

---

## 5. 新机器资源建议

| 项目 | 800 基因、10 网络 × 500 细胞、每重复 |
|---|---|
| CPU | 1 核/任务（**建网实际上单核**，不要给一个任务分配多核） |
| 内存 | 约 2 GB/任务 |
| 单任务耗时 | 1.5–2 小时（含 10 分钟敲除） |

| 机器规格 | 20 个任务（2 条件 × 10 重复）墙钟时间 |
|---|---|
| 8 核 / 16 GB | 约 5–6 小时（并行 6–7） |
| 16 核 / 32 GB | 约 3 小时（并行 12–16） |
| 32 核 / 64 GB | 约 2 小时（并行 20，一轮跑完） |

**并行数请按内存定**：`并行数 × 2 GB < 可用内存`。本机之所以跑不动，就是可用内存只剩 4.2 GB。

若资源仍然紧张，两个降规模方案（需在论文中注明）：

| 方案 | 命令 | 代价 |
|---|---|---|
| 基因数降到 400 | `./run_stability_linux.sh 8 10 400 800` | 与主分析规模（800）不一致，作为补充分析 |
| 重复数降到 5 | `./run_stability_linux.sh 8 5 800 800` | ρ 的估计只有 10 对，置信区间较宽 |

---

## 6. 哪些不需要重跑

包里已包含以下结果，**直接可用**（这也是打包脚本默认带上它们的原因）：

| 分析 | 结果文件 | 重建成本 |
|---|---|---|
| 主分析（两条件网络 + 面板 + 全基因敲除） | `outputs/results/vko_ms_microglia/` | 每网络 2 小时 |
| 技术对照（同条件半样本） | `outputs/results/technical_null/` | 每网络 2 小时 |
| 参数敏感性（4 档） | `outputs/results/param_sensitivity/` | 每档 21 分钟 |
| 跨数据集复现（10 个发表网络） | `results/tables/Table_S3_…csv` | 下载 2 分钟 |
| 四张图 | `results/figures/` | 秒级（需上游结果） |

---

## 7. 跑完后的收尾

1. 把 `Table_S10`、`Table_S11` 拷回主仓库的 `results/tables/`
2. 按检查点结果修改 `docs/09_methods_and_limitations_draft_EN.md`：
   - 四个检查点全过 → **删除 Limitations 第 5 条**，在 Methods 中补一段稳定性描述
   - C2（ρ ≥ 0.6）不通过 → 把"输出不可复现"从 Limitations 提升为 Results 的一个小节
3. `CITATION.cff` 里的 `<user>/<repo>` 与作者信息填好，推送到 GitHub 并归档到 Zenodo 获取 DOI

---

## 8. 当前机器（搬迁前）的状态

| 任务 | 状态 |
|---|---|
| 深度匹配敏感性分析 | 运行中（Control 50% / MS 60%，预计 1–2 小时完成）。可直接中止，也可等它写完再打包——写完后 `outputs/results/vko_ms_depthmatched/` 会有完整结果 |
| 稳定性复现 | 启动过一次，因启动器脚本的 PowerShell 数组解包 bug 中断；**bug 已修复**，并已用单个重复手工验证 worker 脚本可正常运行 |
| 已完成并可打包 | 主分析、技术对照、参数敏感性、跨数据集复现、CellOracle 对照、四张图、全部文档 |

**建议的搬迁顺序**：

1. 等深度匹配那轮写完（或直接中止）
2. 运行 `make_transfer_bundle.ps1 -Full`（如果新机器不便重下数据就用全量）
3. 压缩、拷贝
4. 在新机器上按第 3、4 节操作
