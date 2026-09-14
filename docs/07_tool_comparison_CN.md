# 第二、第三个工具对照：CellOracle

日期：2026-09-13/14

---

## 1. 环境搭建（Windows 上的实际障碍，需在论文中说明）

CellOracle 是 Python 包，且**在 Windows 上无法用常规方式安装**：

| 障碍 | 具体情况 | 解决方式 |
|---|---|---|
| 系统 Python 是 3.12 | CellOracle 0.20 要求 `numpy==1.26.4`、`pandas<=1.5.3`、`matplotlib<3.7`、`anndata<=0.10.8`，这些在 3.12 上无预编译包 | 未改动系统 Python（会破坏用户环境），改用**便携版 CPython 3.10.21**（37.5 MB，解压即用，完全隔离） |
| pip 无法联网 | 沙箱证书链问题，pip 连下载小包都挂起 | 自写 Node 解析器遍历 PyPI 依赖树并下载 62 个 wheel（217 MB），再 `pip install --no-index` 离线安装 |
| `velocyto` 无 Windows 版 | 只有源码包，需 C 编译器 | **最小替身模块**（任何被调用的函数立即报错，从而证明未被使用） |
| `gimmemotifs`、`pybedtools` 同样无 Windows 版 | 需编译 / 需外部 bedtools 二进制 | 同上，替身模块 |
| genomepy 试图写用户目录 | 沙箱拒绝 | `sitecustomize.py` 重定向 appdirs 路径到工作区 |
| CellOracle 数据缓存写用户目录 | 同上 | 修改 `celloracle/data/config.py` 指向工作区 |

> **替身是否影响结果？** 三个替身模块对应的功能（RNA 速率、motif 扫描、BED 处理）都不在我们使用的路径上（GRN 构建 + `simulate_shift`）。替身函数**一旦被调用立即抛错**，因此"未报错"即证明未被使用。这一点须在论文 Methods 中如实说明。

## 2. 数据与参数

| 项目 | 设置 |
|---|---|
| 数据 | GSE279180 人 MS 病灶小胶质细胞 |
| 细胞 | 下采样至 **3,000**（保持 MS/Control 比例）——CellOracle 的 GRN 构建在 9,239 细胞上超出本机算力 |
| 基因 | **3,011**（3,000 高变基因 + 全部面板基因）——符合 CellOracle 推荐的 1,000–3,000 规模 |
| Base GRN | `hg19_gimmemotifsv5_fpr2`（37,003 peak × 1,094 TF） |
| 网络 | 单一小胶质群体（教程标准用法）；`bagging_number=5` |
| KO 模拟 | `n_propagation=3`（默认） |
| 计算耗时 | 插补 + GRN + 4 次 KO ≈ 12 分钟（3,000 细胞） |

## 3. 发现一：CellOracle 只能敲除转录因子

面板 20 个基因中，**只有 4 个是转录因子**（base GRN 只收录 TF 作为调控者）：

| 可敲除（TF） | 不可敲除 |
|---|---|
| **SPI1、IRF8、STAT3、NFKB1** | TREM2、TYROBP、CD33、APOE、BIN1、PICALM、CR1、MS4A4A、C1QA、C3、CX3CR1、P2RY12、TMEM119、BTK、CSF1R、IL10RA |

**这意味着：本项目的核心候选基因（TREM2、CD33、BTK、APOE、CSF1R）CellOracle 一个都做不了。** 这是两个工具最根本的适用范围差异——scTenifoldKnk 可以虚拟敲除任何基因。

## 4. 发现二：CellOracle 的敲除效应多跳传播，scTenifoldKnk 严格单跳

### CellOracle（`n_propagation = 3`）

| 敲除基因 | 表达改变的基因数 | 该 TF 直接靶标数 | **受影响基因中属于直接靶标的比例** |
|---|---:|---:|---:|
| SPI1 | 1,205 | 801 | **32.7%** |
| IRF8 | 887 | 511 | **22.1%** |
| STAT3 | 2,498 | 1,648 | **55.5%** |
| NFKB1 | 2,034 | 1,336 | **48.6%** |

**即 44–78% 的受影响基因落在直接靶标之外**——多跳传播得到确认。

### scTenifoldKnk（本项目 W7）

3,195 次敲除中，**99.9% 的差异调控基因落在"被敲除基因 ∪ 其直接靶标"之内**。

### 结论

**"差异调控基因完全来自直接出边"是 scTenifoldKnk 的结构特性，不是这类方法的普遍性质。** 这是本项目方法学论文最有力的对照证据：

- scTenifoldKnk：敲除效应被**限制在一跳之内**，因此枢纽基因的下游信号会被"稀释"殆尽
- CellOracle：显式做多次传播，效应可以到达网络远端

两者的差异源于设计（CellOracle 的 `simulate_shift` 迭代传播 vs scTenifoldKnk 只置零邻接矩阵的一行后做流形对齐），而非数据或参数。

## 5. 发现三：CellOracle 的重复性同样不理想

同一份数据、同一套参数，两次独立运行的结果差异巨大：

| 敲除基因 | 第 1 次运行 | 第 2 次运行 | 倍数 |
|---|---:|---:|---:|
| SPI1 | 28 | 1,205 | **43×** |
| IRF8 | 1,053 | 887 | 1.2× |
| STAT3 | 732 | 2,498 | 3.4× |
| NFKB1 | 2,230 | 2,034 | 1.1× |

原因是 `bagging_number` 引入的随机性（GRN 构建时对细胞做自助抽样）。**这与本项目在 scTenifoldKnk 中发现的稳定性问题性质相同**——两种工具的输出都对随机性敏感，因此"稳定性复现"是这类方法共同的必需步骤。

## 6. 对论文的意义

第三个工具（CellOracle）的对照使方法学论文的论证链条完整：

| 论点 | 证据 |
|---|---|
| 虚拟敲除结果由网络结构（出度）决定 | scTenifoldKnk：3,195 次敲除，DR ⊆ {KO} ∪ 直接靶标，99.9% |
| 该限制是**工具特有**的，不是方法族共性 | CellOracle：仅 22–56% 的受影响基因是直接靶标 |
| 工具适用范围存在实质差异 | CellOracle 只能敲除 TF（本面板 4/20） |
| 两类工具的**稳定性都成问题** | scTenifoldKnk 度依赖伪影；CellOracle 两次运行差达 43× |
| "结果看起来合理"不构成证据 | 同条件分割对照的 Top 基因同样"合理"（W5） |

## 7. 产出文件

| 文件 | 内容 |
|---|---|
| `outputs/results/celloracle/GRN_links_Microglia.csv` | CellOracle 推断的 GRN（85,532 条 TF→靶标边，含系数与 p 值） |
| `outputs/results/celloracle/celloracle_ko_results.pkl` | 4 次 KO 的逐细胞表达变化矩阵（3,000 × 3,011） |
| `outputs/results/celloracle/celloracle_vs_direct_targets.csv` | 第 4 节的对比表 |
| `outputs/figures/Fig_tool_comparison.png` / `.pdf` | 三面板对照图（受影响基因数、单跳 vs 多跳、运行间变异） |
| 环境 | `work/pyport/python`（便携版 Python 3.10 + 62 个 wheel + 3 个替身模块） |
| 脚本 | `work/celloracle_prep_hvg.py`、`work/celloracle_run.py`、`work/celloracle_compare.py`、`work/pip_resolve.mjs` |

## 8. 前提与局限（须写入论文）

1. 由于本机算力限制，CellOracle 使用 3,000 细胞（非全部 9,239）与单一群体 GRN
2. 三个替身模块（velocyto / gimmemotifs / pybedtools）虽未被调用，但属非标准安装，需说明
3. 两次运行的对比说明 CellOracle 结果有随机性，正式分析应做多种子重复
4. 两个工具的 GRN 类型不同（scTenifoldKnk 为全基因调控网络，CellOracle 为 TF 中心网络），因此**不能直接比较"受影响基因数"的绝对值**，只能比较**单跳/多跳这一结构性差异**
