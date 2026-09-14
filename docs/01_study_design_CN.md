# scTenifoldKnk × 神经免疫疾病：可行性判断与启动方案

> 适用场景：以常规 scRNA-seq / snRNA-seq 数据为基础，在**不重建动物模型、不做真实基因敲除**的前提下，对候选基因做"虚拟敲除"筛选，用于基因功能预测、调控网络解析与治疗靶点优先级排序。
> 版本核实日期：2026-09-13（scTenifoldKnk CRAN v1.1，2026-09-02 发布）

---

## 0. 一句话结论

**可以用，而且是神经免疫领域性价比很高的一步，但它的正确定位是"假设生成器 + 优先级排序器"，不是因果证据。**

它回答的问题只有一个：*在某个细胞类型的调控网络里，把基因 X 的调控输出抹掉，网络里哪些基因/通路会失稳？*
它不回答：X 敲低后细胞表型怎么变、X 在体内的因果作用、X 与其它细胞的通讯如何改变。

神经免疫研究恰好适合它，因为：神经免疫的核心细胞（小胶质细胞、边界相关巨噬细胞、CD4⁺/CD8⁺ T 细胞、B 细胞、星形胶质细胞）都有大量公开的人源单细胞数据，而**人的神经免疫细胞几乎不可能做真实基因敲除**——这正是虚拟敲除的用武之地。

---

## 1. 工具现状（已核实）

| 项目 | 状态 |
|---|---|
| CRAN 正式版 | **scTenifoldKnk v1.1**，2026-09-02 发布，`install.packages("scTenifoldKnk")` |
| 源码仓库 | github.com/cailab-tamu/scTenifoldKnk（286 stars，最近提交 2026-09-02，仍在活跃维护） |
| 许可证 | GPL (≥2) |
| 关键依赖 | `scTenifoldNet (>= 1.4)`、Matrix、MASS、cli、enrichR、igraph、reshape2；可选 `locfdr`（经验零分布） |
| 原始文献 | Osorio D, Zhong Y, …, Cai JJ. *Patterns* 2022;3(3):100434（被引 ≈187） |
| 使用广度 | Europe PMC 中提及该工具的文献 ≈256 篇，其中神经/精神类已覆盖卒中、帕金森病、阿尔茨海默病、癫痫、抑郁症、唐氏综合征、孤独症 |
| 其它语言实现 | Python: scTenifoldpy；MATLAB: scGEAToolbox |
| 姊妹工具 | **scTenifoldXct**（*Cell Systems* 2023，预测细胞–细胞互作）、**scTenifoldNet**（两条件网络比较） |

**一句话现状**：这是一个已进 CRAN、仍在维护、被大量"计算+少量实验验证"型论文采用的成熟工具；同时也已经出现对其可靠性的系统性评测（见第 4 节）。

---

## 2. 它到底模拟了什么（决定你怎么写论文）

### 2.1 技术路径

```
原始 counts（仅需 1 个条件下的细胞）
   ↓ scQC          细胞/基因质控
   ↓ CPM 归一化
   ↓ makeNetworks  对细胞多次子采样 → 主成分回归构建 n 个网络（默认 10 个 × 每次 500 细胞）
   ↓ tensorDecomposition  CP 张量分解去噪得到共识网络
   ↓ strictDirection 强制边的方向性
   ↓ 【虚拟敲除】把目标基因的所有出边（outdegree）置 0
   ↓ manifoldAlignment   野生型网络 vs 敲除网络做非线性流形对齐
   ↓ dRegulation    每个基因在两条件下的欧氏距离 → Box-Cox → Z → χ²(1) 检验 → BH FDR
输出：每个基因的 distance / Z / FC / p.value / p.adj（即"虚拟敲除扰动基因"）
```

### 2.2 它能给出的东西

- 一张**扰动谱**：敲除 X 后，网络中被显著扰动的基因排名（可直接喂给 GSEA/GO/KEGG）
- 一个**节点重要性度量**：X 是不是网络枢纽（hub）
- 一张**以 X 为中心的子网络图**（`plotKO()`）
- 多基因/多条件下的**扰动谱可比较性**：可做聚类、PCA、相关性分析

### 2.3 它明确**不**模拟的东西（论文里必须写清楚）

| 不模拟 | 后果 | 补救方式 |
|---|---|---|
| 剂量效应 | 原生只有"敲除"（出边置 0），没有原生"敲低"参数 | 自行做部分敲低扩展（见第 8 节），标明为未验证扩展 |
| 蛋白层 / 翻译后修饰 / 代谢 | 只在转录调控层面 | 结合蛋白组、功能实验 |
| 配体–受体跨细胞通讯 | 网络局限在单一细胞群内部 | 用 **scTenifoldXct** 或 CellChat/CellPhoneDB 补 |
| 细胞命运转变动力学 | 无时间维度，本质是稳态相关 | 用 RNA velocity / 谱系追踪补 |
| 脱靶与代偿 | 拓扑视角，不建模补偿回路 | 实验对照（scramble / 非靶向 sgRNA） |
| 真实表达量变化 | **虚拟敲除改变的是调控边，不是表达值** | 措辞上永远写 "in silico knockout / predicted perturbation"，不写 "X knockdown decreased Y" |

---

## 3. 神经免疫领域的现成先例（可直接对标）

### 3.1 人 MS 病灶中的真实应用（最有参考价值）

**Maggi P, …, Absinta M. *eBioMedicine* 2023;94:104701**（PMID 37437310，开放获取）

- 数据：已发表的人 MS **慢性活动病灶边缘**免疫细胞 snRNA-seq（GRN 含 6266 个基因）
- 做法：用 scTenifoldNet 模拟"敲除"稀疏淋巴细胞亚群（MS4A1⁺CD20 B 细胞、浆母细胞、T 细胞亚型），并用 **scTenifoldKnk 虚拟敲除 BTK**（BTK 在小胶质和 B 系细胞均表达）
- 预测：B 细胞耗竭会影响小胶质细胞中与**铁/血红素代谢、缺氧、抗原提呈**相关的基因
- 验证：72 例 MS 患者、约 2 年 MRI 随访（46 例抗 CD20 治疗 vs 26 例未治疗）→ **抗 CD20 并未使顺磁环（PRL）消退**
- 意义：这是"计算预测被真实世界数据部分否定"的范例写法。**审稿人喜欢这种诚实**，而不是硬把预测说成真理

### 3.2 工具作者自带的神经相关案例（代码可直接复用）

GitHub 仓库 `inst/manuscript/` 下现成的分析：

| 目录 | 内容 | 对你的价值 |
|---|---|---|
| `TREM2/` | GSE130626（cuprizone 小鼠小胶质细胞，WT vs *Trem2* KO）；**只用 WT 细胞建网做虚拟敲除**，再与真实 KO 的差异基因对照 | 「用真实 KO 数据评估虚拟敲除」的标准模板 |
| `MECP2/` | Rett 综合征相关 *MECP2* 虚拟敲除 | 神经发育–免疫交叉先例 |
| `MICROGLIA/` | 微胶质细胞数据整合脚本 | 数据准备参考 |
| `STABILITY/` | 10 次 500 细胞子采样复现、单/双/三基因敲除（Trem2、Trem2+Apoe、Trem2+Apoe+Lpl） | **稳定性分析的官方做法** |
| `NEGCONTROL/` | 以 Cftr/Akap7/Zranb1/Krcc1/Mta1/Rps8 为敲除对象做阴性对照 | **阴性对照设计模板** |

### 3.3 已发表的其它神经/精神类应用

- 缺血性卒中单核转录组 + 铁死亡虚拟敲除（*Heliyon* 2024）
- 帕金森病空间代谢组 + 单细胞虚拟敲除筛选（*Metabolites* 2026）
- 阿尔茨海默病蛋白稳态失衡与星形胶质细胞 KLF15–LRRTM4 轴（2026 预印本）
- 癫痫星形胶质细胞色氨酸代谢（*Front Neurosci* 2026）
- 抑郁症岩藻糖基化 / 星形胶质细胞致病模式（*Front Neurosci* 2026）
- 唐氏综合征中孕期皮层细胞类型特异特征（*Nat Commun* 2025）
- 孤独症 IGF1R 多组学（*Front Genet* 2024）

**结论：神经免疫 + 虚拟敲除不是空白区，但"跨细胞类型 / 跨疾病条件 / 带实验验证"的组合仍然稀缺——这是你的切入口。**

---

## 4. 必须知道的负面证据（2026 预印本，务必读）

**Wu S, …, Ge W. "Method Choice, Not Biology, Determines In Silico Perturbation Results: A Systematic Evaluation of Eight Methods Across Four Datasets." 2026 预印本**（DOI 10.64898/2026.08.11.744106，尚未同行评审）

要点：

1. 8 种方法（覆盖 6 类数学框架）× 4 个数据集的系统评测中，**8 种方法里有 6 种（含 VAE 类与张量分解类）检测不到 TF→通路 的方向性信号**；相对表现较好的是 CellOracle 与 DDIM。
2. **DDIM 与 scTenifoldKnk 的基因排序显著负相关（ρ = −0.811, p = 0.027）——换一个方法就能翻转生物学结论。**
3. CellOracle 预测的扰动方向与 CRISPRi Perturb-seq 实验方向一致率仅 40.9%（与随机无差异），说明**稳态相关性 ≠ 因果扰动**。
4. 作者给出的最低数据建议：**≥500 个细胞、≥1000 个高变基因**。
5. 失败模式诊断：图方法的主要问题是"图非特异性"（富集度 0.84×，即接近背景）。

**你的应对策略（直接写进 Methods 和 Limitations）**

- 至少用**两种原理不同的方法**交叉（如 scTenifoldKnk + CellOracle 或 scRank），报告一致性/不一致性，而不是只报一种
- 必做**稳定性分析**（10 次子采样）与**阴性对照基因**
- 对外宣称的结果**必须落到可实验验证的 3–5 个基因**，不要给 200 个基因的清单
- 结论措辞：*predicted / prioritized / hypothesis-generating*，而非 *demonstrated*

---

## 5. 六个递进式选题（筛选 → 验证 → 机制 → 体内）

以模板「基因 X × 神经免疫疾病 Y × 细胞类型 Z」代入（示例：`TREM2 × 多发性硬化 × 病灶边缘小胶质细胞`）。

**题目 1｜建立 Z 细胞类型特异的单细胞调控网络：数据整合、注释校正与网络可复现性**
→ 用公开 scRNA/snRNA-seq 数据（多供体），IUCT/Seurat/Harmony 整合，人工标志基因 + 参考映射双重注释，输出一个"细胞数、基因数、供体来源、注释置信度"全部可追溯的建网输入矩阵。

**题目 2｜基因 X 及候选面板的虚拟敲除筛查：扰动谱驱动的靶点优先级排序**
→ scTenifoldKnk 面板虚拟敲除（20 个候选 + 4 个阴性对照），输出扰动谱、Z 值排序、核心被扰动模块，并给出排序的稳定性（10 次子采样 Spearman ρ）与阴性对照零分布。

**题目 3｜疾病与对照网络中的扰动谱比较：识别疾病特异性易感节点**
→ 对照与疾病条件分别建网，同一面板各跑一次虚拟敲除，比较扰动谱差异（相似度矩阵、PCA、置换检验），找出"只在疾病网络中失稳"的基因——这是相对单条件分析最容易做出新意的模块。

**题目 4｜从细胞自主到细胞间通讯：虚拟敲除后的小胶质–淋巴/星形胶质交互重塑**
→ scTenifoldXct 预测敲除后配体–受体互作变化，叠加 CellChat/NicheNet，得到"基因 X → 小胶质 → T 细胞/星形胶质"的跨细胞假说链。

**题目 5｜计算预测的实验验证：iPSC 来源小胶质细胞 CRISPRi 敲低与多重功能读出**
→ iPSC 分化小胶质细胞（或原代人小胶质）CRISPRi/siRNA 敲低，读出细胞因子分泌（IL-6/TNF-α/IL-10）、吞噬、ROS、铁代谢、MHC-II 提呈、迁移；与虚拟敲除预测的被扰动通路做一致性检验。

**题目 6｜体内因果与转化价值：模型动物中的条件性敲低与表型救援**
→ EAE（MS）、cuprizone 脱髓鞘、5xFAD/APP 模型中的小胶质条件性敲低（Cx3cr1-CreER 等），评估疾病评分、脱髓鞘/髓鞘再生、免疫浸润；结合人类遗传学证据（GWAS/eQTL）判断转化价值。

> 逻辑链：**网络可复现（1）→ 计算筛选（2）→ 疾病特异性（3）→ 机制扩展（4）→ 实验验证（5）→ 体内因果（6）**。
> 若只做计算类论文（无湿实验），建议止步于 1–4，并把 4 作为讨论与局限，避免声称机制已证实。

---

## 6. 推荐主线与三个具体 Aims

**推荐组合：题目 2 为骨干 + 题目 3 作为核心增量 + 题目 4 作为机制扩展 + 题目 5 的 3–5 个基因验证。**

**核心科学问题（一句话）**
> 在疾病 Y 的 Z 细胞类型调控网络中，候选基因集合 S 的虚拟敲除会扰动哪些核心调控模块，其中哪些基因的扰动最具疾病特异性、稳定且可被实验验证？

**Aim 1｜构建可复现的细胞类型特异 scGRN**
多供体数据整合 → 细胞类型/亚群注释 → 对照与疾病两套建网输入矩阵 → 网络稳健性报告（细胞数、基因数、供体贡献、参数敏感性）。

**Aim 2｜面板虚拟敲除与优先级排序**
20 个候选基因 + 4 个阴性对照 → `transcriptomeWide` 一次建网、多次扰动 → 扰动谱统计 → 稳定性（10× 子采样）→ 阴性对照零分布 → 疾病 vs 对照扰动谱差异 → 排序出 Top 5。

**Aim 3｜三重验证**
① 计算复现：独立数据集/独立供体分层；② 外部队列对照：与已发表真实 KO（如 GSE130626 的 *Trem2* KO）、Perturb-seq/CRISPRi 公共数据比较；③ 湿实验：iPSC-小胶质 CRISPRi 敲低 + 功能读出。

---

## 7. 数据准备（决定成败的一步）

### 7.1 输入要求

| 项目 | 要求 | 说明 |
|---|---|---|
| 输入格式 | **原始 counts**，基因 symbol 为行、细胞为列 | `qc = TRUE` 时**不要预先归一化** |
| 细胞数 | 每个细胞类型 **≥500**（建议 ≥1000） | 官方子采样默认每次 500 细胞；低于 500 时网络不稳定 |
| 基因数 | 1000–5000（HVG 或表达过滤后基因） | 基因数直接影响运行时间（见 7.4） |
| 细胞类型 | **必须按细胞类型分开建网** | 混合细胞类型会得到被异质性污染的"平均网络" |
| 条件 | 对照与疾病各一套网络 | 这是设计上的核心增量 |
| 物种 | 人/鼠均可，注意基因 symbol 体系一致（hsa2mmu 需转换） | 换算会丢失部分基因 |

### 7.2 数据来源建议

- **GEO**：检索 `(disease) AND (single-cell OR single-nucleus) AND (microglia OR CSF OR lesion OR brain)`
- **CELLxGENE Census**：按细胞类型（microglial cell、T cell…）和疾病直接切片，元数据规范
- **AD Knowledge Portal / Synapse**：SEA-AD、ROSMAP 等大规模 AD 小胶质数据
- **Human Cell Atlas / 单细胞脑图谱**：如 Science 2024 的 388 人脑调控网络资源（适合做外部网络对照）
- **已发表 MS 病灶 snRNA-seq**（如 Maggi/Absinta 系列）与 **EAE 小鼠队列**（做跨物种一致性）

### 7.3 数据清单模板（建议先填满再动手）

| GSE/数据集 | 物种 | 组织 | 平台 | 细胞数 | 目标细胞类型数 | 条件 | 供体数 | 注释质量 | 是否用于复现 |
|---|---|---|---|---|---|---|---|---|---|
| | | | | | | | | | |

### 7.4 运行时间预估（官方基准，单核）

| 细胞数 | 基因数 | 一次完整单基因虚拟敲除 |
|---:|---:|---:|
| 300 | 1,000 | ≈3.5 min |
| 1,000 | 1,000 | ≈4.3 min |
| 1,000 | 5,000 | ≈2 h 52 min |
| 5,000 | 5,000 | ≈3 h 9 min |
| 7,500 | 7,500 | ≈10 h 16 min |

**关键点**：耗时几乎全部在**建网**这一步。因此：

- 20 个基因的筛选请用 `transcriptomeWide = TRUE`：**只建一次网，然后逐个基因扰动**，额外成本主要是流形对齐
- 10 次稳定性复现 ≈ 10 次建网 ≈ 数十小时量级 → 用 `nCores` 并行，或把面板降到 5–8 个基因、基因数控制在 2000–3000

---

## 8. 运行流程与关键参数

配套脚本：`outputs/scripts/run_vko_panel.R`（可直接改配置运行）

### Step 0｜环境

```r
install.packages("scTenifoldKnk")          # v1.1 来自 CRAN
install.packages(c("locfdr", "fgsea", "igraph", "enrichR", "Matrix"))
# 数据处理：Seurat / Bioconductor 按需
```

### Step 1｜数据准备与注释

读取 counts → Ensembl 转 symbol → 读取 metadata（barcode, celltype, condition, donor）→ **亚聚类确认目标细胞类型纯度**（至少检查标志基因与双细胞评分）。

### Step 2｜质控（官方 `scQC` 行为）

- `qc_minLibSize = 1000`（细胞最低文库）
- `qc_removeOutlierCells = TRUE`（1.58×IQR/√n 离群细胞剔除）
- `qc_maxMTratio = 0.1`（线粒体比例；注意基因名需为 `MT-` 前缀）
- `qc_minPCT = 0.05`（基因至少在 5% 细胞中表达）
- **额外过滤**：`^Rp[0-9]|^Rpl|^Rps|^Mt-`（核糖体与线粒体基因，官方脚本做法）

### Step 3｜建网参数（默认值与调参建议）

| 参数 | 默认 | 建议 |
|---|---|---|
| `nc_nNet` | 10 | ≥10；网络数越多越稳，代价线性上升 |
| `nc_nCells` | 500 | 取 min(500, 细胞数×0.8)；细胞少时降低并报告 |
| `nc_nComp` | 3 | 3–10，做敏感性分析 |
| `nc_q` | 0.9 | 0.9–0.95，控制网络稀疏度 |
| `nc_lambda` | 0 | 0–0.5，控制方向性强度；官方对 TREM2 做过 0/0.25/0.5/0.75/1 的敏感性分析 |
| `td_K` | 3 | 3–5（CP 张量秩） |
| `nc_priorNetwork` | NULL | **可选的加分项**：传入已知 TF–靶基因先验（如 SCENIC regulon、ChIP/ATAC 证据）约束网络 |
| `dr_empiricalNull` | FALSE | 设为 **TRUE**（需 locfdr）用经验零分布，对"全局大面积扰动"的假阳性更稳健，建议两种都跑并报告 |

### Step 4｜虚拟敲除的两种模式

```r
# 模式 A：单基因（返回完整 diffRegulation 表）
res <- scTenifoldKnk(countMatrix = X, gKO = "TREM2")

# 模式 B：面板 / 转录组范围（建网一次，多次扰动）
tw <- scTenifoldKnk(countMatrix = X, transcriptomeWide = TRUE,
                    gKO = c(panel_genes, neg_controls))
tw$perturbationDistances   # 行 = 被敲基因，列 = 全部基因
```

### Step 5｜转录组范围模式的统计量（必须自己算）

模式 B 只返回距离矩阵，需要按官方 `dRegulation()` 的同一公式换算：

```r
d  <- tw$perturbationDistances["TREM2", ]
FC <- d^2 / mean(d^2)
p  <- pchisq(FC, df = 1, lower.tail = FALSE)
p.adj <- p.adjust(p, method = "fdr")
# Z 来自距离的 Box-Cox 变换后再标准化（脚本中已实现，与 scTenifoldKnk::dRegulation 一致）
```

> 一致性检查：建议先对 1–2 个基因分别跑模式 A 与模式 B，确认两种路径得到的 Z/p 排序高度一致，再批量跑面板。

### Step 6｜稳定性分析（官方做法，必做）

复刻 `inst/manuscript/STABILITY/`：对 10 次不同随机种子、每次 500 细胞的子采样重复全流程 → 收集每个基因的 Z 向量 → 分位数归一化 → 计算两两 **Spearman ρ** → 报告均值±标准差与分布图。
作者的 TREM2 结果以相关矩阵 + 山脊图（ridges）呈现，可直接照搬图型。

### Step 7｜阴性对照（必做）

选 4 个**在该细胞类型中几乎不表达或生物学无关**的基因（官方示例用了 Zranb1、Krcc1、Mta1、Rps8 等）一起虚拟敲除：

- 期望：它们的扰动谱彼此不相关、且显著基因数很少
- 若阴性对照也给出上百个显著基因 → 说明网络/参数或细胞数有问题，**先修数据，别继续解释结果**

### Step 8｜功能解读

- 对 Z 值排序做 **fgsea**（KEGG/Reactome/GO）或 `enrichR`
- 用 `plotKO(res, gKO = "TREM2")` 画以敲除基因为中心的子网络（可叠加富集类别饼图）
- 与疾病 GWAS 风险基因、eQTL、已知通路做一致性检验

### Step 9｜方法交叉（应对第 4 节的风险）

至少再上一个原理不同的方法（CellOracle / scRank / SCENIC+ / Dictys / GEARS），比较：
Top-N 基因重叠率（Jaccard）、Z 值相关（Spearman）、通路富集一致性。
**方法间不一致本身就是一个值得报告的结果**，不要藏起来。

### Step 10｜结果产出清单

扰动距离矩阵 · 每个基因的 DR 表 · 面板汇总排序表 · 稳定性相关矩阵 · 阴性对照零分布 · 富集表 · 子网络图 · 疾病 vs 对照扰动谱对比图 · 参数敏感性热图 · 运行参数与 sessionInfo。

---

## 9. "敲低（knockdown）"而不是"敲除（knockout）"怎么办

**原生工具没有部分敲低参数**——它只有"把出边全部置零"。若你的课题确实要模拟敲低（例如 siRNA/CRISPRi 而非 KO），可以做一个显式标注为**扩展**的剂量分析：

```r
WT <- res$tensorNetworks$WT
KO <- WT
KO["TREM2", ] <- KO["TREM2", ] * (1 - kd)     # kd = 0.25 / 0.5 / 0.8
MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
DR <- scTenifoldKnk::dRegulation(MA)           # 支持 empiricalNull = TRUE
```

建议同时跑 `kd = 0, 0.25, 0.5, 0.8, 1.0`，绘制**剂量–网络响应曲线**：

- 单调递增 → 说明扰动效应稳健，是加分项
- 非单调或跳变 → 提示网络对度值敏感，需要更保守的结论

论文中必须写明："partial perturbation was implemented by scaling the gene's outdegree edges, which is an extension of the original binary knockout scheme and was not part of the validated workflow."

---

## 10. 验证的分层结构（写论文时按此分层，避免过度声称）

| 层级 | 做了什么 | 能证明 | **不能**证明 |
|---|---|---|---|
| 网络构建 | 多细胞子采样 + 张量分解 + 稳定性复现 | 网络在给定数据与参数下可复现 | 网络是真实生物学网络 |
| 虚拟敲除 | 扰动谱 + 负控 + 剂量曲线 | 该基因在**计算网络**中是重要节点 | 真实细胞中该基因的功能与因果 |
| 富集/网络 | 通路与子网络 | 功能关联假说 | 机制因果 |
| 外部数据对照 | 与真实 KO/Perturb-seq 一致 | 预测具有一定的外部效度 | 在你的细胞体系里成立 |
| 细胞实验 | iPSC-小胶质 CRISPRi + 功能读出 | 该基因在**该系统**中调控所测表型 | 体内因果、临床疗效 |
| 动物实验 | EAE/cuprizone/AD 模型条件性敲低 | 体内因果与疾病表型贡献 | 人类临床转化 |

---

## 11. 审稿人风险清单（提前写好 Limitations）

1. **方法可靠性**：2026 预印本显示同类方法间结论可相互矛盾（ρ = −0.811）→ 双方法 + 稳定性 + 实验验证；引用该预印本并说明其未经同行评审。
2. **网络推断的敏感性**：结果依赖细胞数/HVG 选择/`nc_q`/`td_K`/`nc_lambda` → 做参数网格敏感性分析并报告，而非只报一组"好看"的参数。
3. **细胞类型异质性**：混合细胞群或注释错误会污染网络 → 亚聚类质控、双细胞检测、标志基因核对。
4. **虚拟敲除 ≠ 真实敲低**：不建模剂量、代偿、蛋白层 → 措辞严格、加剂量曲线、湿实验验证。
5. **无法捕捉细胞间通讯**：单细胞群内部网络 → scTenifoldXct / CellChat 互补。
6. **泛化性不足**：单数据集单供体 → 多队列/多供体分层复现。
7. **富集结果过度解读**：图方法本身富集特异性偏低（0.84×）→ 报告零分布对照（作者脚本中的 shuffle 检验做法）。
8. **阴性对照缺失**：审稿人常问"你怎么知道这不是网络噪声？" → 第 7 节的负控零分布是直接答案。

---

## 12. 12 周启动时间表

| 周次 | 任务 | 交付物 |
|---|---|---|
| W0 | 环境安装（R + scTenifoldKnk v1.1 + locfdr/fgsea）、跑通官方示例、跑通脚本 `quick_test` 模式 | 可运行的本地环境 |
| W1 | 候选基因清单定稿（20 个候选 + 4 个阴性对照）、数据检索与清单填写、细胞类型注释 | 数据清单表 + 注释后的对象 |
| W2 | 目标细胞类型亚聚类确认、QC、建网输入矩阵（对照/疾病各一套） | 两套 `countMatrix` + QC 图 |
| W3 | 小规模试跑（300 细胞 × 1000 基因，1 个基因）→ 模式 A 与模式 B 一致性核对 | 试跑报告 |
| W4 | 全量建网 + 24 基因面板虚拟敲除（transcriptomeWide） | 距离矩阵 + 面板 DR 表 |
| W5 | 稳定性分析（10× 子采样）+ 阴性对照 + 参数敏感性 | 稳定性 ρ 矩阵、负控零分布 |
| W6 | 疾病 vs 对照扰动谱比较；富集；子网络图 | 全部主图 |
| W7 | 第二方法交叉（CellOracle/scRank）与一致性分析 | 方法一致性报告 |
| W8 | 与外部真实 KO/Perturb-seq 数据对照；Top 5 定稿 | 验证优先级表 |
| W9 | 实验设计定稿（iPSC-小胶质 CRISPRi、读出指标、样本量） | 实验方案 |
| W10–12 | 启动湿实验（分化/敲低/表型读出）；同步撰写方法学与结果 | 初稿 + 实验数据 |

---

## 13. 第一周就能做的 6 件事

1. `install.packages("scTenifoldKnk")`，跑官方 README 里的负二项模拟示例（2000 细胞 × 100 基因），确认输出结构与 `plotKO()` 可用
2. 打开作者仓库的 `inst/manuscript/TREM2/Code/TREM2_DataProcessing.R`，逐行理解"用 WT 建网 + 与真实 KO 对照"的完整链路
3. 复刻 `inst/manuscript/STABILITY/stabilityAnalysis.R` 的稳定性框架（换成你自己的基因）
4. 复刻 `inst/manuscript/NEGCONTROL/negControls.R` 的阴性对照框架
5. 确定你的 20 候选基因来源（GWAS 风险基因 / 差异表达 hub 基因 / 药物靶点 / 文献候选），并为每个基因写一句"为什么值得预测"
6. 选定 1 个公开数据集，跑通"注释 → QC → 300 细胞 × 1000 基因小规模建网"，用真实数据感受一次运行时间

---

## 14. 资源、文献与引用

**工具与代码**

- CRAN: https://cran.r-project.org/package=scTenifoldKnk （v1.1, 2026-09-02）
- GitHub: https://github.com/cailab-tamu/scTenifoldKnk （含 TREM2 / MECP2 / MICROGLIA / STABILITY / NEGCONTROL 全套分析脚本）
- scTenifoldNet（两条件网络比较）: https://github.com/cailab-tamu/scTenifoldNet
- scTenifoldXct（细胞–细胞互作）: https://github.com/cailab-tamu/scTenifoldXct
- Python 实现 scTenifoldpy / MATLAB 实现 scGEAToolbox

**必引文献**

1. Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Tian Y, Chapkin R, Huang JZ, Cai JJ. **scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation.** *Patterns*. 2022;3(3):100434. doi:10.1016/j.patter.2022.100434
2. Osorio D, et al. **scTenifoldNet: A machine learning workflow for constructing and comparing transcriptome-wide gene regulatory networks from single-cell data.** *Patterns*. 2020;1(9):100139. doi:10.1016/j.patter.2020.100139
3. Maggi P, …, Absinta M. **B cell depletion therapy does not resolve chronic active multiple sclerosis lesions.** *eBioMedicine*. 2023;94:104701. doi:10.1016/j.ebiom.2023.104701（神经免疫领域的 scTenifoldKnk 应用范例）
4. Wu S, Hu G, Yang Z, Wang Z, Cai J, Mao J, Ge W. **Method Choice, Not Biology, Determines In Silico Perturbation Results: A Systematic Evaluation of Eight Methods Across Four Datasets.** 2026 预印本. doi:10.64898/2026.08.11.744106（**方法学风险，必读**）
5. Kamimoto K, et al. **Dissecting cell identity via network inference and in silico gene perturbation.** *Nature*. 2023;614:742–751. doi:10.1038/s41586-022-05688-9（CellOracle，方法交叉的首选）
6. Qiu R, Zhao MM. **A confound-diagnostic toolkit for in silico perturbation with single-cell foundation models.** 2026 预印本. doi:10.64898/2026.08.04.732812（混淆因素诊断，补充讨论用）

---

## 15. 配套脚本

`outputs/scripts/run_vko_panel.R` 提供：

- 数据读取（10x MTX / RDS / CSV）与元数据筛选
- 官方同款 QC 与核糖体/线粒体基因过滤
- 面板虚拟敲除（`transcriptomeWide = TRUE`，一次建网多次扰动）
- 与官方 `dRegulation()` 完全一致的 Z / FC / p / FDR 换算
- 多次子采样稳定性分析（Spearman ρ + 排名一致性）
- 阴性对照零分布
- fgsea 富集与 `plotKO` 子网络输出
- 结果表、图、参数与 `sessionInfo` 全量落盘

**重要提示**：脚本中的参数块需按你的数据修改；脚本在无 R 环境的机器上编写，**首次运行请先用 `quick_test = TRUE` 小规模验证**，再跑全量。
