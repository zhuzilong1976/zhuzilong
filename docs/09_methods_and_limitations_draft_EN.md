# Methods and Limitations — draft

> 说明（中文）：这是基于本项目 W1–W9 全部实测结果写成的论文 Methods + Limitations 初稿，正文为英文以便直接投稿使用。
> 所有数字都可追溯到 `outputs/W1`–`W9` 文档与 `outputs/results/` 下的结果文件；带 **[待补]** 标记处需要在投稿前补充。

---

## Methods

### Data

**Primary dataset.** We used single-nucleus RNA-seq of microglia from human subcortical multiple sclerosis (MS) lesions (GEO accession **GSE279180**, "Cell type mapping reveals tissue niches and interactions in subcortical multiple sclerosis lesions"). The study provides cell-type–split AnnData files; we used the microglia subset (`GSE279180_ctype_MG.h5ad`), which contains **23,627 genes × 9,239 microglia** as raw integer counts. The metadata distinguish condition (Control, *n* = 2,005 cells; MS, *n* = 7,234), lesion type (chronic active CA, *n* = 4,410; chronic inactive CI, *n* = 2,824; control tissue Ctrl, *n* = 2,005), 13 donors, and eight microglial subtypes. Median library size was 927 counts in Control and 3,564 in MS cells.

**Choice of contrast.** We initially considered a within-disease contrast (chronic active vs chronic inactive lesions), which is better matched for sequencing depth (median library size 3,888 vs 3,144; cell-number ratio stable at 1.55–1.62 across filtering thresholds). We rejected this design because **the two lesion classes are derived from completely non-overlapping donors** (CA: MS197, MS229, MS377, MS411; CI: MS497, MS549; with MS549 alone contributing 83.1% of CI cells), so any lesion-type effect would be inseparable from donor effects. We therefore used a case-control contrast (MS vs Control) and report the depth imbalance explicitly (see Limitations).

**Query gene panel.** We defined a panel of 20 candidate genes related to myeloid risk and neuroimmune regulation (TREM2, TYROBP, CD33, APOE, BIN1, PICALM, CR1, MS4A4A, SPI1, IRF8, STAT3, NFKB1, C1QA, C3, CX3CR1, P2RY12, TMEM119, BTK, CSF1R, IL10RA) plus four negative-control genes (ACTA2, TSHR, PECAM1, YKT6). All 20 candidates were expressed above the detection threshold in the combined microglia data; expression was strongly MS-enriched for most (e.g. APOE 72.5% vs 23.8% of cells; TYROBP +22.0 percentage points; C1QA +23.5), except C3, which was lower in MS (69.8% vs 77.8%).

**Independent replication datasets.** To test the generality of our findings we used the network objects published by the developers of scTenifoldKnk in their own manuscript folder (`inst/manuscript/*/Results/*.RData` of the `cailab-tamu/scTenifoldKnk` repository), comprising **10 independent networks** spanning human and mouse, five tissues (liver, lung, muscle, intestine, brain), network sizes of 2,591–9,970 genes, and seven knocked-out genes (Hnf4a+Hnf4g, Nkx2-1, Dmd, Malat1, Ahr, Mecp2 ×2, Cftr ×2, Akap7). These were generated with earlier versions of the software and processed with the developers' own pipelines; we used them as-is.

### Network construction and virtual knockout

Networks were built with **scTenifoldKnk 1.1** (CRAN, released 2026-09-02) on R 4.6.1, following the package pipeline: quality control → CPM normalisation → construction of *n* principal-component-regression networks from subsampled cells → CANDECOMP/PARAFAC (CP) tensor decomposition → strict directionality → transposition, after which the adjacency-matrix rows correspond to regulators.

Parameters for the primary analysis were the package defaults except where stated: `nc_nNet = 10`, `nc_nCells = 500` (bootstrap sampling **with replacement**, as implemented in `scTenifoldNet::makeNetworks`), `nc_nComp = 3`, `nc_q = 0.90`, `nc_lambda = 0`, `td_K = 3`, `ma_nDim = 2`. Quality control used `scQC` with `minLibSize = 500`, `removeOutlierCells = TRUE`, `maxMTratio = 0.1`, and ribosomal/mitochondrial genes (`^Rp[0-9]|^Rpl|^Rps|^Mt-`) removed. Because `scQC` applies its gene-expression filter per condition, we replaced it with a **shared gene set** computed on all target cells (expression in >5% of the 9,239 microglia) so that both condition-specific networks contain identical genes; the 800 most variable genes of this shared set were retained, with all panel genes forced to be included. This left 800 genes and, after cell-level QC, **1,347 Control cells** and **6,927 MS cells**.

A virtual knockout was performed for **every gene in the network** (not only the panel): the gene's entire row (outgoing edges) was set to zero, a non-linear manifold alignment was computed between the wild-type and knocked-out networks (`scTenifoldNet::manifoldAlignment`, *d* = 2), and differential regulation was tested (`scTenifoldKnk::dRegulation`). We computed statistics both with the theoretical χ²(1) null and, where noted, with Efron's empirical null (`locfdr`). For each knockout *g* we defined the **differentially regulated set** DR(*g*) = {genes with Benjamini–Hochberg adjusted *p* < 0.05}, and the **outdegree** of *g* as the number of non-zero entries in its row of the wild-type adjacency matrix.

### Containment analysis

For each knockout *g* we tested whether

> DR(*g*) ⊆ {*g*} ∪ {direct out-neighbours of *g* in the wild-type network}

and whether the weaker bound |DR(*g*)| ≤ outdegree(*g*) + 1 held. We applied this test to **all genes** of four networks built from the MS-microglia data (Control network, MS network, and two control split-half networks; 800 + 800 + 797 + 798 = **3,195 knockouts**) and to the 10 published replication networks. In the replication set the knocked-out gene was identified directly from the data as the row(s) containing non-zero edges in `WT` that are all-zero in `KO`, avoiding reliance on file naming.

### Technical control (split-half null)

To estimate the false-positive rate of the two-condition comparison, we randomly split the Control cells into two halves (446 and 450 cells after quality control), built **two independent networks from the same biological condition** with identical parameters, and performed the same manifold alignment and differential-regulation analysis. This comparison contains no biological difference; its significant-gene rate is therefore the technical null.

### Condition comparison

We compared the Control and MS networks directly with `scTenifoldNet::manifoldAlignment` followed by `scTenifoldKnk::dRegulation` — the two-condition mode for which the scTenifold framework was originally designed — and contrasted the resulting significant-gene rate with the split-half null (Fisher's exact test).

### Empirical null for the knockdown analysis

To place the panel genes' effects on an empirical scale rather than an arbitrary threshold, we computed a specificity score for every gene. Let *d*(*k*, *j*) be the distance of target *j* under knockout *k*. We standardised each target across all 800 knockouts of a network,

> *z*(*k*, *j*) = ( *d*(*k*, *j*) − mean_k *d*(·, *j*) ) / sd_k *d*(·, *j*) ,

and summarised each knockout by max_z = max_{j≠k} *z*(*k*, *j*) and by *n*_{z>2} = #{*j* : *z*(*k*, *j*) > 2}. Because max_z is computed over ~800 targets, its null distribution is not centred at zero; we therefore used the **776 non-panel genes as an empirical null** and report empirical *p*-values with Benjamini–Hochberg correction, plus a version restricted to genes with outdegree within ±20% of the query gene.

### Parameter sensitivity

To test whether the observed behaviour is a property of the package defaults (`nc_lambda = 0`), we rebuilt 400-gene Control networks under four settings — `nc_lambda` = 0, 0.5, 1.0 (the latter enforcing strictly directional edges) at `nc_q` = 0.90, and `nc_lambda` = 0 at `nc_q` = 0.95 — each with 10 networks × 400 cells, and repeated the full all-gene knockout analysis.

### Stability analysis

To separate the structural properties of the virtual knockout from properties of one particular cell sample, we repeated the entire single-condition analysis on **10 independent cell subsamples per condition**. For replicate *r* of a condition we drew 800 cells at random from that condition (seed 20260913 + *r*), applied the same quality-control pipeline as in the primary analysis (`scQC`, `minLibSize = 500`, `removeOutlierCells = TRUE`, `maxMTratio = 0.1`, ribosomal and mitochondrial genes removed), kept the same shared gene set (the 800 most variable genes of the combined microglia data, panel genes forced in), rebuilt 10 networks from 500 cells each (`q = 0.90`, `K = 3`, `nc_nComp = 3`) and knocked out **every** gene of the network (797–800 knockouts per replicate; 20 replicates in total).

Four checkpoints were specified before the runs: **(C1)** containment (DR ⊆ {knockout} ∪ direct out-neighbours) in ≥ 99% of knockouts per replicate, excluding knockouts with outdegree 0, whose knockout is a no-op; **(C2)** median pairwise Spearman ρ of the per-gene max_z ranking between replicates of the same condition ≥ 0.6; **(C3)** zero panel genes at FDR < 0.05 under the empirical null (non-panel knockouts, Benjamini–Hochberg) in every replicate; **(C4)** no gene at FDR < 0.05 in the two-condition comparison in any replicate. For C4 the comparison was repeated with Control replicate *r* paired with MS replicate *r* (manifold alignment, *d* = 2, followed by `dRegulation`), so that both networks of a pair were built from an independent cell sample.

### Comparison with a second method (CellOracle)

We compared our findings with **CellOracle 0.20.0**, a widely used in-silico knockout method of different design (it constructs a TF-centric GRN from a motif-based base GRN and propagates the perturbation iteratively). Because CellOracle cannot be installed on Windows through standard channels, we constructed an isolated environment: a portable CPython 3.10.21 distribution, 62 dependency wheels resolved and downloaded independently, and minimal placeholder modules for three packages with no Windows distribution (`velocyto`, `gimmemotifs`, `pybedtools`). The placeholders raise an exception upon any call; no exception occurred during the analyses reported here, i.e. the code paths used (GRN construction and `simulate_shift`) do not depend on them.

We used the human promoter base GRN (`hg19_gimmemotifsv5_fpr2`; 37,003 peaks × 1,094 TFs) and the microglia data reduced to the standard CellOracle scale: 3,000 cells (subsampled preserving the MS/Control ratio) and 3,011 genes (3,000 highly variable genes plus all panel genes). Networks were built with `bagging_number = 5` on a single microglial GRN unit, and knockouts simulated with `n_propagation = 3` (default). For each knockout we compared the set of genes with non-zero expression change (`delta_X`) with the TF's direct targets in the inferred GRN, and performed the whole analysis twice to gauge run-to-run variability.

### Statistics and software

All *p*-values are two-sided. Multiple testing was controlled with the Benjamini–Hochberg procedure. Hypergeometric tests were used for pathway over-representation with the 800 or 3,011 network genes as background. Enrichment libraries were KEGG 2021, Reactome 2022, GO Biological Process 2021 and MSigDB Hallmark 2020. Comparisons of significant-gene rates used Fisher's exact test. Analyses were performed in R 4.6.1 (`Matrix`, `scTenifoldKnk` 1.1, `scTenifoldNet` 1.4, `locfdr`, `fgsea`, `igraph`) and Python 3.10.21 (`scanpy` 1.10.4, `anndata` 0.10.8, `CellOracle` 0.20.0). Scripts are listed under Data and code availability.

### Data and code availability

All data are public (GEO: GSE279180; the published replication networks are distributed with the scTenifoldKnk source repository). The analysis scripts, including the stability replicates and the checkpoint evaluation (`code/03_controls/28_`, `29_`, `code/07_supplementary/67_stability_checkpoints_C3C4.R`), accompany this repository. The derived data — complete distance matrices for all knockouts, the 20 stability replicate objects, the split-half control networks and the CellOracle results — are archived separately because of their size (Zenodo DOI **[待补: DOI]**).

---

## Results

### Stability of the knockout and comparison outputs across independent cell subsamples

All 20 replicates (10 per condition) completed; after quality control they contained 343–369 Control and 669–699 MS cells, with network densities of 84.6–88.5% and 77.1–82.1%, respectively (Table S10, Figure S2A). Per replicate, 1–13 genes had outdegree 0, so their knockout is a no-op.

**Containment is reproducible (C1).** Excluding knockouts with outdegree 0, containment held for **100% of knockouts in all 20 replicates** (the strict value including outdegree-0 knockouts was 99.25–99.87%). The containment rule is therefore a property of the network construction rather than of the particular cells sampled.

**The ranking of perturbation magnitude is not (C2).** The Spearman correlation between the max_z rankings of two replicates of the same condition had a **median of 0.185** over all 90 replicate pairs (range −0.466 to 0.818), well below the pre-specified threshold of 0.6 (Figure S2B). Which targets appear most strongly perturbed by a given knockout therefore depends strongly on the cell subsample.

**Panel significance under the empirical null is reproducible (C3).** In every one of the 20 replicates, **no** candidate panel gene reached FDR < 0.05 against the empirical null built from the non-panel knockouts; the smallest candidate FDR across all replicates was 0.062 (Figure S2C). The conclusion that the panel is indistinguishable from random network genes is stable.

**The two-condition comparison is not (C4).** Pairing Control replicate *r* with MS replicate *r*, **6 of 10 pairs** produced at least one gene at FDR < 0.05 (1–2 genes per pair; smallest FDR 0.0069, recurring genes CD74, ADGRB3, CHST11, MAN2A1, TTC7A; Figure S2D), whereas the primary analysis reported none. The number of genes at *p* < 0.05 was 41–54 per pair and the correlation between the two networks of a pair was 0.22–0.29. The absence of significant genes in the single primary comparison is therefore not a stable property of the data.

Together, C1 and C3 show that the *qualitative* findings (locality of the perturbation; panel genes indistinguishable from random network genes) survive resampling, whereas C2 and C4 show that magnitudes — both the ranking of perturbed targets and the exact number of genes crossing an FDR threshold in the condition comparison — do not. Biological reading of individual target genes, or of a "no significant genes" result, should therefore rest on averages across replicates rather than on a single network.

**Caveat on the size of the replicates.** Each replicate is built from 800 sampled cells (343–369 Control and 669–699 MS cells after quality control), i.e. less than half to a quarter of the corresponding primary network (1,347 and 6,927 cells). Part of the variability measured by C2 and C4 therefore reflects the smaller sample size rather than a property of the primary networks; the asymmetry between the stable structural results (C1, C3) and the unstable magnitude results (C2, C4) is nevertheless reproduced in both conditions and at both replicate sizes.

---

## Limitations

1. **Scale.** Because of local compute constraints, the primary networks contain 800 genes. The developers' own analyses use 2,000–10,000 genes, and the API advises 1,000–3,000 genes; we therefore cannot exclude that larger networks would show a qualitatively different pattern of differential regulation. The parameter-sensitivity experiments were run at 400 genes for the same reason, although the property tested is structural and the replication set spans 2,591–9,970 genes.

2. **Single dataset for the application.** The MS-microglia application uses one study. Its case-control contrast is affected by a **substantial depth imbalance** (median library size 1,700 in Control vs 3,517 in MS after quality control; depth ratio 2.07). We mitigated this in the primary analysis by lowering the cell-level library-size threshold to 500 and by reporting per-sample composition. **Depth-matched sensitivity analysis.** We additionally rebuilt both networks restricting cells to a common depth window, which we selected by scanning candidate windows for the joint objective of retaining sufficient cells and equalising the depth distributions: the window 1,500–5,000 counts reduced the median-depth ratio from 2.07 to ~1.20 and the Kolmogorov–Smirnov statistic between the two depth distributions by 12 orders of magnitude (KS *p* from 2.1 × 10⁻¹²⁶ to 1.2 × 10⁻¹⁴), while retaining 752 Control and 5,762 MS cells after quality control. Full depth matching is not achievable because the two groups' depth distributions have different shapes.

3. **Donor composition.** MS cells are dominated by two donors (MS549 32.5%, MS377 29.1%), and the Control and MS groups share no donors — an inherent property of case-control designs, but one that limits the inference that observed differences are disease-driven.

4. **The containment rule says nothing about biological validity.** Our results characterise where differentially regulated genes come from (the knocked-out gene's direct out-neighbours) and how their number is governed (outdegree, via weight dilution). They do **not** show that the implicated targets are functionally important. In the MS-microglia data the candidate panel was indistinguishable from random network genes under an empirical null, and a same-condition split-half comparison produced as many "significant" genes as the disease comparison. In addition, the number of genes crossing FDR < 0.05 in the two-condition comparison varies between replicate pairs (0–2 genes), so the absence of significant genes in the primary comparison should not be over-interpreted.

5. **Stability.** The qualitative conclusions replicate across independent cell subsamples: with 10 replicates per condition (800 genes, 10 networks × 500 cells each), containment held for 100% of non-degenerate knockouts in all 20 replicates, and no candidate panel gene reached FDR < 0.05 against the empirical null in any replicate. Magnitude-based statements do not replicate: the median Spearman correlation of the max_z ranking between replicates was 0.185 (range −0.466 to 0.818), and 6 of 10 replicate-paired condition comparisons produced 1–2 genes at FDR < 0.05 while the primary comparison produced none. We therefore report the containment rule and the panel result as robust, and treat all per-gene effect sizes, rankings and "not significant" statements about the condition comparison as sample-dependent.

6. **The second-method comparison is not like-for-like.** CellOracle builds a TF-centric GRN from a motif-based base GRN and requires the knocked-out gene to be a transcription factor (only 4 of our 20 panel genes are), whereas scTenifoldKnk can knock out any gene. CellOracle was run on 3,000 cells and a single GRN unit. The two methods can therefore be compared with respect to the *locality* of the perturbation (one-hop versus multi-hop), but not with respect to absolute numbers of affected genes.

7. **Software environment.** The CellOracle results were produced with three dependency placeholders (see Methods). Although the placeholders are never invoked in the reported analyses, this is a non-standard installation and should be reproduced in a native environment (e.g. Linux) before publication.

8. **Mixed species and annotation.** The replication networks mix human and mouse data and were processed with heterogeneous pipelines and software versions; gene symbols were used as provided, without harmonisation. This heterogeneity strengthens the generality of the containment finding but prevents quantitative cross-dataset comparison of effect sizes.

9. **Negative controls.** The negative-control genes available for this dataset (ACTA2, TSHR, PECAM1, YKT6) are expressed at moderate levels but have low network degrees, which itself affects how their knockouts behave. Degree-matched negative controls would be preferable in future work.

10. **Empirical-null dependence of degenerate knockouts.** For genes with zero outdegree the knockout is a no-op and the resulting statistics are unstable: the same analysis yields 0 significant genes under a theoretical χ² null in one network and 188 under an empirical null in another. Such genes should be excluded rather than reported, and our numbers for them should not be interpreted.

---

## 中文小结（供你决定下一步）

**Methods 已经具备完整可写的材料**：数据、设计取舍（为什么放弃 CA vs CI）、网络参数、包含关系的定义、三类对照（技术对照 / 经验零分布 / 参数敏感性 / 稳定性复现）、跨数据集复现、第二工具对照、统计与软件，全部有实测数字支撑。

**2026-09-14 更新：稳定性复现已在新机器上完成**，Methods 与 Results 两节已按实测数字写好（见上文 "Stability analysis" 与 "Stability of the knockout and comparison outputs"）。结果：

- C1 通过：20/20 重复的包含关系合规率 100%（剔除出度 0 的敲除）；
- **C2 不通过**：max_z 排名的重复间 Spearman ρ 中位 0.185（阈值 0.6）；
- C3 通过：20/20 重复中没有任何面板基因经验零分布 FDR<0.05（最小 0.062）；
- **C4 不通过**：10 对重复的条件比较中有 6 对出现 1–2 个 FDR<0.05 的基因。

Limitations 第 4、5 条已按上述结果改写，不再有 [待补] 标记。

**仍未完成的事项**：

1. 公开代码仓库与 DOI（`Data and code availability` 里的 [待补]）；
2. C4 的脚本已归档为 `code/07_supplementary/67_stability_checkpoints_C3C4.R`，并已写进 Data and code availability；
3. 软件版本差异一句话说明：稳定性复现用的是 R 4.6.1 + scTenifoldKnk 1.1 + scTenifoldNet 1.4（与主分析一致），仅 `fgsea` 版本不同（1.38.0 vs 1.34.2），且未重跑富集分析；
4. 深度比的两个基准（"2.07 → 1.20" 是 QC 后细胞；"3.84 → 1.20" 是全部细胞）需统一表述。
