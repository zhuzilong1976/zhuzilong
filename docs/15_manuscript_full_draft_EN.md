# Manuscript draft (full) — target journal: Briefings in Bioinformatics

> **状态（2026-09-14）**：完整初稿，所有数字均来自本仓库的结果文件（最后一节给出"主张 → 证据文件"对照表）。
> 仍需作者补充的位置一律标 `[TO BE COMPLETED]`；无法核验的文献标 `[CITATION NEEDED]`，未凭空生成任何参考文献。
> 引用的 9 篇文献已全部通过 CrossRef 逐条核验（标题/期刊/年份/卷页/DOI）。
> 文章类型：**Original Article**（2026-09-14 确定）。

---

## Title

**Containment, not biology: the differentially regulated gene set of single-cell virtual knockout is the knockout's own out-neighbourhood**

**Running title:** What single-cell virtual knockout actually measures

## Authors and affiliations

Zhu Zilong^1,*^

^1^ Department of Neurology, Tianjin Huanhu Hospital, Tianjin, China

\* Correspondence: zhuzilong1976@gmail.com ｜ ORCID: 0000-0002-6955-0903

---

## Abstract

Virtual knockout on single-cell gene regulatory networks is increasingly used to nominate genes and drug targets, and the resulting "differentially regulated" (DR) gene list is usually read as the downstream consequence of losing a gene. We show that this list is, instead, a structural property of the network. Across **3,195 knockouts in four human microglia networks**, the DR set was contained in the knocked-out gene together with its **direct out-neighbours** in **99.87%** of cases; all four exceptions were knockouts with **outdegree 0**, for which the knockout is a numerical no-op yet which returned up to 15 significant genes. Under the package's default χ² null, **96% of knockouts returned exactly one significant gene — the knocked-out gene itself**; under the empirical null the list expands to **outdegree + 1** while the knockout's weight is concentrated on its targets (198 significant genes at outdegree 197) and collapses to a single gene once that weight is diluted (outdegree ≥ 584). In the multiple sclerosis network the **two negative-control genes produced the two largest DR sets** (188 and 146 genes; p~adj~ down to 10⁻²¹⁶) because one of them had outdegree 0. The rule replicated in **10 independently published networks** (2,591–9,970 genes, human and mouse, five tissues; 100% containment) and in the developers' own published *Trem2* result (128/128 significant genes inside the direct target set). It is **tool-specific**: with CellOracle, which propagates perturbations iteratively, only 22–56% of affected genes were direct targets. The practical consequence is that a disease-versus-control network comparison must be read against a same-condition null: the MS-versus-control rate (7.25% of genes at *p* < 0.05) did not exceed a split-half null built from the same condition (5.65%; Fisher *p* = 0.22), no gene passed FDR control in either, and the null's top genes were as biologically plausible as the disease comparison's; the candidate panel was indistinguishable from random network genes. Across 10 independent cell subsamples per condition, containment and the panel result replicated exactly, whereas magnitude-based outputs did not (median ranking correlation 0.185; 6 of 10 replicate-paired comparisons returned 1–2 genes at FDR < 0.05). We give four reporting requirements that make such analyses interpretable.

**Keywords:** in-silico knockout; gene regulatory network; scTenifoldKnk; CellOracle; reproducibility; outdegree; single-cell transcriptomics

## Key points

* The "differentially regulated" output of single-cell virtual knockout is **contained in the knocked-out gene and its direct out-neighbours** (99.87% of 3,195 knockouts; 100% in 10 published networks and in the developers' own *Trem2* result).
* The **number** of genes in that list is governed by the knockout's **outdegree** — not by biology. In the MS network, two negative-control genes topped the panel (188 and 146 genes) because one had outdegree 0.
* The default χ² null returns essentially one gene (the knockout itself) in 96% of knockouts; the empirical null returns outdegree + 1. The chosen null must be reported.
* A disease-versus-control comparison must be accompanied by a **same-condition null**: the significant rate (7.25%) did not exceed the split-half null (5.65%; Fisher *p* = 0.22), and the null's top genes were equally plausible.
* Qualitative results (containment; panel indistinguishable from random genes) replicate across cell subsamples; magnitude-based results (rankings; which genes cross FDR) do not.

---

## 1. Introduction

Single-cell gene regulatory networks (scGRNs) are now used to simulate the loss of a gene in silico. In the most widely used implementation, `scTenifoldKnk`, a network is inferred from single-cell transcriptomes, the row of the target gene is zeroed (its outgoing edges removed), and the wild-type and knocked-out networks are compared by non-linear manifold alignment; genes whose distance changes beyond a significance threshold are reported as "differentially regulated" (DR) (Osorio *et al.*, 2022). The output is a ranked gene list, and the field's working assumption is that this list describes what happens downstream of losing the gene (Maggi *et al.*, 2023). A second family of tools takes a different route: `CellOracle` builds a transcription-factor-centric network from a motif-based base GRN and propagates a perturbation iteratively (Kamimoto *et al.*, 2023).

Two developments make the interpretation of these outputs urgent. First, systematic benchmarking has shown that in-silico perturbation results are dominated by method choice rather than biology: eight methods evaluated across four datasets produced largely method-specific gene lists (Wu *et al.*, 2026, preprint). Second, the confounders of these pipelines are beginning to be catalogued, including the dependence of perturbation effects on network structure (Qiu and Zhao, 2026, preprint).

Here we characterise, at the level of the mechanism, what the DR set of `scTenifoldKnk` actually contains. Using human multiple sclerosis (MS) lesion microglia as a worked example, and then replicating in published networks from five tissues and two species, we show that:

1. the DR set is **contained** in the knocked-out gene plus its direct out-neighbours — it is not a downstream readout;
2. the **size** of the DR set is set by the knocked-out gene's **outdegree**, through dilution of each target's incoming weight, so that biologically irrelevant genes with low or zero outdegree can produce the largest lists;
3. consequently, the "plausibility" of a virtual-knockout gene list, and the outcome of a disease-versus-control network comparison, cannot be interpreted without a same-condition null;
4. the containment property replicates across parameter settings, network sizes, species and datasets, but magnitude-based outputs do not survive cell subsampling.

We close with four reporting requirements that we argue should accompany any virtual-knockout analysis.

---

## 2. Results

### 2.1 The differentially regulated set is contained in the knockout's own out-neighbourhood

We inferred networks from the microglia subset of a single-nucleus atlas of subcortical MS lesions (GSE279180; 23,627 genes × 9,239 cells; associated publication *Nature Neuroscience* 2024, doi:10.1038/s41593-024-01796-z). Four networks were built with identical parameters (10 networks × 500 cells, *q* = 0.90, *K* = 3): the Control network (800 genes; 1,347 cells after QC), the MS network (800 genes; 6,927 cells), and two split-halves of the Control cells (797 genes/446 cells and 798 genes/450 cells). Every gene of each network was knocked out, giving **3,195 knockouts**; for each knockout *g* we tested whether its DR set was a subset of {*g*} ∪ {direct out-neighbours of *g*} (Figure 1A, Table S1).

Containment held for **3,191 of 3,195 knockouts (99.87%)**, and it held for **100% (all 3,191) of knockouts with outdegree > 0**. All four violations were knockouts with **outdegree 0** — genes whose entire row in the inferred network is zero, so that "knocking them out" changes nothing: `TSHR` in the Control network (15 significant genes), `ACTA2` and `FAF1` in the MS network (7 each) and `TSHR` in split-half A (9). The bound |DR(*g*)| ≤ outdegree(*g*) + 1 followed directly.

An independent check inside the developers' own data gave the same picture: in their published *Trem2* mouse microglia result, *Trem2* has outdegree 2,744 and 128 significant genes, **all 128 of which were direct targets of *Trem2*** (100%), even though the naive prediction |DR| = outdegree + 1 would have implied thousands.

### 2.2 Under the default null the DR set is essentially the knockout itself

With the package's default theoretical χ²(1) null, the DR set was minute: **3,064 of 3,195 knockouts (95.9%) returned exactly one significant gene**, the median |DR| was 1, and the maximum over all knockouts was 15. The single significant gene was, in the great majority of cases, the knocked-out gene itself. The proportion of direct targets that reached significance was 0.92 for the six knockouts with outdegree 1–10 and **0.000 for every knockout with outdegree > 10** (Figure 1B).

### 2.3 Under the empirical null the DR size tracks outdegree

Because an empirical null gives a different calibration, we repeated the panel analysis with Efron's empirical null (`locfdr`; Tables S7, S8), which makes the outdegree dependence visible. In the Control network, `BTK` (outdegree 197) returned **198** significant genes and `CD33` (outdegree 146) returned **147** — outdegree + 1 in both cases — whereas hub genes with outdegree 682–788 (`TREM2`, `APOE`, `BIN1`, `C1QA`, `C3`, `CR1`, …) returned **1**. The same pattern appeared in the MS network (`YKT6`, outdegree 145 → 146; `CD33`, outdegree 75 → 76).

The interpretation is weight dilution: a target is called significant in proportion to the fraction of its total incoming weight contributed by the knocked-out gene (Figure 1B). When the knockout contributes most of a target's incoming weight (low outdegree), essentially all of its targets are called significant, so |DR| ≈ outdegree + 1; when the knockout is one of hundreds of regulators (high outdegree), its weight is diluted below the significance threshold in every target, and only the knockout itself remains in the list.

### 2.4 Consequence: biologically irrelevant knockouts can top the list

The MS network contained a natural experiment. Two of the four negative-control genes (chosen *a priori* as biologically unrelated to MS microglia) behaved as follows: **`ACTA2` produced the largest DR set in the whole panel — 188 genes, with adjusted *p*-values as low as 1.3 × 10⁻²¹⁶ — and `YKT6` produced the second largest, 146 genes**, because `ACTA2` had outdegree 0 (a numerical no-op) and `YKT6` outdegree 145 in that network. In the Control network the same two genes returned 0, where `ACTA2` had outdegree 500. Two of the four negative controls thus ranked **first and second** by list size in one network and last in the other, with no change in biology. Under the theoretical χ² null, by contrast, the outdegree-0 knockouts returned 7–15 genes — an instability of two orders of magnitude for the same gene, driven purely by the choice of null model.

### 2.5 The rule replicates across published networks and parameter settings

We applied the containment test to **10 network objects published by the developers of `scTenifoldKnk`** (2,591–9,970 genes; human and mouse; liver, lung, muscle, intestine and brain; knocked-out genes *Hnf4a*+*Hnf4g*, *Nkx2-1*, *Dmd*, *Malat1*, *Ahr*, *Mecp2* ×2, *Cftr* ×2, and the negative control *Akap7*; outdegrees from 1 to 3,632). **Containment was 100% in all 10 datasets** (`n_DR` between 2 and 377; Table S3, Figure 2A, B). The bound also held in every case, and the developers' *Trem2* result provided a further independent point with outdegree 2,744 and |DR| = 128 — far below the naive outdegree + 1 expectation, exactly as predicted by weight dilution.

Rebuilding 400-gene Control networks under four parameter settings (*nc_lambda* = 0, 0.5 and 1 at *nc_q* = 0.90; *nc_lambda* = 0 at *nc_q* = 0.95; 403 knockouts each) changed network density from 56.7% to 86.8% and median outdegree from 257 to 373, but not the rule: compliance was 99.0% (4 outdegree-0 violations) under the first three settings and 96.0% (16 outdegree-0 violations) at *q* = 0.95; excluding outdegree-0 knockouts, **compliance was 100% in all four settings** (Table S2, Figure S1).

### 2.6 The property is tool-specific, not general to network-based perturbation

We repeated the analysis with `CellOracle` (0.20.0) on the same microglia data (3,000 cells, 3,011 genes, motif-based base GRN `hg19_gimmemotifsv5_fpr2`, *n*~propagation~ = 3, bagging 5), and compared each knockout's affected-gene set with the transcription factor's direct targets in the inferred GRN (Table S9, Figure 3). Because the base GRN is transcription-factor-centric, only **4 of our 20 panel genes were eligible for knockout**. For those, the fraction of affected genes that were direct targets ranged from **22.1% (*IRF8*) to 55.5% (*STAT3*)** — that is, **44–78% of the affected genes lay outside the direct target set**, consistent with iterative multi-hop propagation. The containment property is therefore a specific signature of the `scTenifoldKnk` implementation (row-zeroing plus one-shot manifold alignment), not of in-silico perturbation in general.

### 2.7 A disease-versus-control comparison without a same-condition null is not interpretable

The framework's original use case is not single-gene knockout but network comparison between conditions. Comparing the Control and MS networks gave **58 of 800 genes at *p* < 0.05 (7.25%)** and **no gene at FDR < 0.05**. To ask whether that excess is signal, we built two networks from randomly split halves of the *same* condition (446 and 450 cells; no biological difference) and repeated the identical analysis: **45 genes (5.65%) at *p* < 0.05**, again none at FDR < 0.05. The difference between the disease comparison and the same-condition null was not significant (**Fisher exact *p* = 0.221**, odds ratio 1.31, 95% CI 0.86–2.00), and the top-100 lists of the two comparisons shared only a Jaccard index of 0.136. Critically, the null's top genes (`SPP1`, `CX3CR1`, `CSF2RA`, `SLC1A3`, `SORL1`) are as plausible as microglial/immune genes as the disease comparison's (`MAN2A1`, `C3`, `ADGRB3`, …). Had we omitted the null, we would have reported a biologically convincing gene list that the pipeline generates from pure noise.

The same conclusion holds for the panel. Under the empirical null, **no candidate panel gene reached FDR < 0.05 in either network** (minimum adjusted *p* 0.525 in Control; Tables S4, S5), and in the condition comparison the panel's mean rank percentile was 48.9% (random genes: 50.1%). The candidate panel was indistinguishable from random network genes.

### 2.8 Depth matching does not change the conclusion

Control and MS cells differ substantially in sequencing depth (median 927 versus 3,564 counts; ratio 3.84). Restricting both conditions to the window 1,500–5,000 counts reduced the median-depth ratio to 1.20 and retained 752 Control and 5,762 MS cells after quality control; the depth-matched networks (798 and 800 genes) were compared with the same pipeline. The two networks were more similar than in the primary analysis (Spearman correlation 0.317 versus 0.196), the number of nominally significant genes remained small (and **0 of 800 at FDR < 0.05**), and the panel's median *p*-value moved from 0.454 to 0.595 (Supplementary Data 1). The absence of a FDR-controlled signal is therefore not an artefact of the depth imbalance.

### 2.9 Which results survive resampling, and which do not

Finally we repeated the entire single-condition analysis on **10 independent cell subsamples per condition** (800 cells sampled per replicate, identical QC and gene selection, 10 networks × 500 cells, every gene knocked out; 20 replicates, 15,971 knockouts in total) (Table S10, Figure S2).

* **Containment replicated exactly**: 100% of non-degenerate knockouts in **all 20 replicates** (1–13 outdegree-0 knockouts per replicate).
* **The panel result replicated exactly**: in **no** replicate did any candidate panel gene reach FDR < 0.05 under the empirical null (smallest candidate adjusted *p* across all replicates: 0.062).
* **Magnitude-based outputs did not replicate**: the median Spearman correlation between the max-*z* rankings of two replicates of the same condition was **0.185** (90 pairs, range −0.466 to 0.818), far below the pre-specified threshold of 0.6.
* **The condition comparison did not replicate**: pairing Control replicate *i* with MS replicate *i* (both networks built from independent cell samples), **6 of 10 pairs returned 1–2 genes at FDR < 0.05** (smallest FDR 0.0069; recurring genes *CD74*, *ADGRB3*, *CHST11*, *MAN2A1*, *TTC7A), whereas the primary comparison returned none; 41–54 genes per pair were significant at *p* < 0.05.

Qualitative statements about the behaviour of the pipeline therefore generalise across cell samples; statements about effect sizes, rankings, and "no significant genes" do not.

---

## 3. Discussion

### 3.1 What the DR set is

Our results identify the mechanism that generates the output of `scTenifoldKnk`. Zeroing the knocked-out gene's row removes its outgoing edges; the subsequent comparison of wild-type and knocked-out networks measures, for every target, how much its position in the manifold alignment changed. A target moves appreciably only if the knocked-out gene supplies a substantial share of its incoming regulatory weight. Two consequences follow. First, no gene outside {*g*} ∪ *N*⁺(*g*) can move — the perturbation is one hop by construction, so the DR set is contained in the knockout's out-neighbourhood (99.87% of 3,195 knockouts; 100% of non-degenerate knockouts). Second, **within** that neighbourhood the number of genes called significant is a monotone function of how concentrated the knockout's weight is: near-complete for low outdegree, negligible for hubs (Figure 1B).

This explains several observations that have been attributed to biology. The developers' *Trem2* knockout returns 128 genes from an outdegree of 2,744; we show those 128 genes are exactly the subset of *Trem2*'s direct targets whose incoming weight is dominated by *Trem2*. In our MS microglia data, a knockout with **outdegree 0** — a numerical no-op — returned 188 "significant" genes with adjusted *p* down to 10⁻²¹⁶, more than any genuine candidate gene, and the panel's two largest lists came from the negative controls. Neither number reflects downstream biology; both are consequences of network structure and of the null model in use.

### 3.2 Relation to recent benchmarking and confound analyses

Systematic benchmarking across eight in-silico perturbation methods concluded that results are determined by method choice rather than by the underlying biology (Wu *et al.*, 2026, preprint). Our work supplies the mechanism for one of the most widely used members of that family: the DR set is a local, weight-filtered property of the inferred network, so any biological interpretation of its content is unsupported unless the network's local structure is accounted for. Complementary work has begun to catalogue confounders of perturbation analyses with foundation models (Qiu and Zhao, 2026, preprint); the outdegree dependence we quantify here is a concrete, reportable instance of such a confounder, and it is measurable from the network alone, before any biological question is asked.

### 3.3 The contrast with CellOracle is structural, not quantitative

With CellOracle, 44–78% of the genes whose expression changed after a knockout were outside the transcription factor's direct target set, because that method propagates the perturbation iteratively. The two tools therefore answer structurally different questions: `scTenifoldKnk` reports the first hop (filtered by weight), CellOracle reports a multi-hop cascade. This implies that comparing the *number* of affected genes between the two tools is meaningless; only structural properties — locality, propagation depth, and the dependence of output size on network topology — are comparable. It also means that "agreement between tools" is not an appropriate validation criterion for a gene list.

### 3.4 Implications for practice: four reporting requirements

**R1. Report the outdegree of every knocked-out gene, and exclude outdegree-0 knockouts.** In our data all containment violations, and both of the largest spuriously large lists, came from genes with outdegree 0; for them the knockout is a no-op. Outdegree is trivially available from the adjacency matrix and should be a required element of any reported knockout.

**R2. Do not rank genes by the number of significant targets.** |DR| is a function of outdegree (outdegree + 1 while weight is concentrated; ≈1 once it is diluted), so ranking by list size ranks network topology, not biology — as our negative controls demonstrate. If a ranking is needed, it must be based on target-level statistics with the outdegree dependence removed.

**R3. Pair every condition comparison with a same-condition null.** Our disease-versus-control comparison produced a 7.25% significant rate that was statistically indistinguishable from the 5.65% produced by comparing two halves of the same condition (Fisher *p* = 0.22), and the null's top genes were equally plausible. A plausible-looking list is therefore not evidence (see also the two-state comparison in Section 2.7). The split-half null costs one extra network and one extra comparison.

**R4. State which null model was used, and report FDR-controlled results.** Under the theoretical χ² null, 96% of knockouts returned a single gene (the knockout itself); under the empirical null, the same knockouts returned outdegree + 1 genes, and degenerate knockouts drifted by two orders of magnitude (7–15 versus 188 genes). Results reported without naming the null, or without FDR control, cannot be reproduced or interpreted.

### 3.5 What this means for the MS application presented here

Applied to human MS lesion microglia, the pipeline produced no FDR-controlled signal between conditions (0 of 800 genes), a candidate panel indistinguishable from random network genes, and a condition-comparison rate that did not differ from a same-condition null. Our study therefore does **not** nominate MS candidate genes; its contribution is to show that this dataset, analysed this way, cannot do so, and to give the diagnostics that make that visible. This is the same conclusion reached by the developers' own published *Trem2* analysis when read structurally, and it is consistent with the negative-control reasoning we recommend.

---

## 4. Limitations

1. **Network scale.** Primary networks contain 800 genes, below the 1,000–3,000 genes advised by the API and below the developers' own 2,000–10,000. The containment property is structural and held in networks up to 9,970 genes, but the equality |DR| ≈ outdegree + 1 was characterised in 800-gene networks; larger networks could shift where the weight-dilution transition occurs.
2. **One dataset for the application.** The MS-microglia application uses a single study with a substantial depth imbalance (median 927 versus 3,564 counts). Depth matching (window 1,500–5,000 counts; ratio 1.20) left the conclusions unchanged, but full matching is impossible because the two depth distributions differ in shape.
3. **Donor composition.** The two conditions share no donors (Control 6 donors, MS 7), and two MS donors contribute most MS cells; any condition effect is therefore inseparable from donor effects — one reason we report no biological claim.
4. **Containment is not validity.** We characterise where the DR genes come from and how many there are; we do not show that the implicated targets are functionally important.
5. **Stability is property-dependent.** Containment and the panel result replicated in all 20 subsample replicates; the ranking of perturbation magnitudes did not (median ρ = 0.185), and 6 of 10 replicate-paired condition comparisons returned 1–2 genes at FDR < 0.05, whereas the primary comparison returned none. The replicates are smaller than the primary networks (343–369 versus 1,347 Control cells), so part of this variability reflects sample size; the asymmetry between structural and magnitude-based results nonetheless appeared in both conditions.
6. **The cross-tool comparison is not like-for-like.** `CellOracle` requires a transcription factor (4 of our 20 panel genes), was run on 3,000 cells and one GRN unit, and its environment was built with three placeholder modules for packages without Windows builds (never invoked). Only the locality of the perturbation is compared, not absolute gene counts.
7. **Replication set heterogeneity.** The 10 published networks mix human and mouse, five tissues and several software versions, with gene symbols used as provided. This strengthens the generality of the containment result but prevents quantitative cross-dataset comparison of effect sizes.
8. **Negative controls.** Our negative-control genes are expressed at moderate levels but have low network degree, and this is precisely why they produce extreme outputs. Degree-matched negative controls would isolate the effect further.
9. **Enrichment analyses were not performed** on the primary networks; the panel's non-enrichment documented here is a null result, not evidence of absence of function.

---

## 5. Conclusions

The "differentially regulated" gene set produced by single-cell virtual knockout is the knocked-out gene's own out-neighbourhood, filtered by how much regulatory weight the knockout contributes to each target. It is contained in that neighbourhood (99.87% of 3,195 knockouts; 100% of non-degenerate knockouts; 100% in 10 published networks), its size is governed by outdegree, and biologically irrelevant knockouts can therefore produce the largest lists. The property is specific to this implementation (CellOracle propagates multi-hop), robust to parameters and depth, and reproducible across cell subsamples for qualitative but not magnitude-based outputs. Interpreting these outputs as downstream biology requires, at minimum, reporting outdegree, avoiding list-size rankings, pairing every condition comparison with a same-condition null, and naming the null model.

---

## 6. Methods

### 6.1 Data and pre-processing

We used the microglia subset of a published single-nucleus atlas of subcortical MS lesions (GEO **GSE279180**; associated publication: *Nature Neuroscience* 2024, doi:10.1038/s41593-024-01796-z), file `GSE279180_ctype_MG.h5ad`: **23,627 genes × 9,239 microglia** as raw counts, with 2,005 Control cells (6 donors) and 7,234 MS cells (7 donors), 13 donors in total and eight annotated microglial subtypes. Median library size was 927 counts (Control) and 3,564 (MS). We rejected a within-disease contrast (chronic active versus chronic inactive lesions), which is better depth-matched, because the two lesion classes are derived from non-overlapping donors, making lesion type inseparable from donor identity.

Quality control followed the package's `scQC` (`minLibSize = 500`, `removeOutlierCells = TRUE`, `maxMTratio = 0.1`), with ribosomal and mitochondrial genes removed, retaining 1,347 Control and 6,927 MS cells. Because `scQC` filters genes per condition, we instead used a **shared gene set** computed across all microglia (expression in > 5% of cells; 7,779 genes), from which the 800 most variable genes were retained with all panel genes forced in, so that both condition-specific networks contain identical genes.

### 6.2 Network construction and virtual knockout

Networks were built with `scTenifoldKnk` 1.1 (CRAN, 2026-09-02) on R 4.6.1, following the package pipeline: CPM normalisation → *n* principal-component-regression networks from subsampled cells → CANDECOMP/PARAFAC tensor decomposition → strict directionality → transposition, so that rows correspond to regulators. Defaults were used except where stated: `nc_nNet = 10`, `nc_nCells = 500` (bootstrap with replacement), `nc_nComp = 3`, `nc_q = 0.90`, `nc_lambda = 0`, `td_K = 3`, `ma_nDim = 2`. A virtual knockout sets the target gene's entire row to zero; the wild-type and knocked-out networks are compared with `scTenifoldNet::manifoldAlignment` (*d* = 2) and `scTenifoldKnk::dRegulation`. Unless noted, statistics use the package's theoretical χ²(1) null; where indicated we used Efron's empirical null (`locfdr`). For each knockout *g*, DR(*g*) = {genes with Benjamini–Hochberg adjusted *p* < 0.05} and outdegree(*g*) = number of non-zero entries in row *g* of the wild-type adjacency matrix (Benjamini and Hochberg, 1995; Efron, 2004).

### 6.3 Containment analysis

For every gene in each network we tested DR(*g*) ⊆ {*g*} ∪ *N*⁺(*g*) and |DR(*g*)| ≤ outdegree(*g*) + 1. We applied the test to all genes of four MS-microglia networks (Control, MS, and two Control split-halves: 800 + 800 + 797 + 798 = **3,195 knockouts**) and to the 10 published replication networks, where the knocked-out gene was identified from the data as the row(s) that are non-zero in `WT` and all-zero in `KO`. Parameter sensitivity used four 400-gene Control networks (403 knockouts each).

### 6.4 Empirical null, specificity and the panel

For each knockout we standardised each target across all knockouts of a network, *z*(*k*, *j*) = (*d*(*k*, *j*) − mean~k~*d*(·, *j*)) / sd~k~*d*(·, *j*), and summarised each knockout by max-*z* and by *n*~z>2~. Because max-*z* is a maximum over ~800 targets, we used the 776 non-panel genes as an empirical null and report empirical *p*-values with Benjamini–Hochberg correction, plus an outdegree-matched variant (±20%) (Tables S4–S8; Efron's empirical null).

### 6.5 Same-condition (split-half) null and condition comparison

Control cells were split at random into two halves (446 and 450 cells after QC), two networks were built with identical parameters, and the same manifold-alignment and differential-regulation analysis was applied. The significant-gene rate was compared with the disease-versus-control comparison by Fisher's exact test. This null contains no biological difference, so its rate estimates the pipeline's false-positive behaviour.

### 6.6 Depth-matched sensitivity

Both conditions were restricted to cells with 1,500–5,000 counts (selected by scanning candidate windows for cell yield and depth equalisation), then processed identically: 798 genes × 752 Control cells and 800 genes × 5,762 MS cells were retained. The comparison of the two depth-matched networks was repeated as in 6.5 and contrasted with the primary networks (`panel_depthmatched_vs_main.csv`).

### 6.7 Stability analysis

The entire single-condition pipeline was repeated on 10 independent cell subsamples per condition (800 cells sampled per replicate with seeds 20260913 + *r*; identical QC, shared gene set, 10 networks × 500 cells, q = 0.90, K = 3; every gene knocked out; 15,971 knockouts in total). Four checkpoints were specified before the runs: **C1** containment ≥ 99% per replicate excluding outdegree-0 knockouts; **C2** median pairwise Spearman ρ of the max-*z* ranking between replicates ≥ 0.6; **C3** zero panel genes at FDR < 0.05 under the empirical null in every replicate; **C4** no gene at FDR < 0.05 in any replicate-paired condition comparison (Control replicate *i* versus MS replicate *i*).

### 6.8 Cross-tool comparison

`CellOracle` 0.20.0 was installed in an isolated CPython 3.10.21 environment (three dependency placeholders, never invoked). We used the human promoter base GRN `hg19_gimmemotifsv5_fpr2` (37,003 peaks × 1,094 transcription factors) and 3,011 genes (3,000 highly variable genes plus the panel), a single GRN unit with `bagging_number = 5`, and knockouts simulated with `n_propagation = 3`. For each transcription factor we compared the set of genes with non-zero expression change with the factor's direct targets in the inferred GRN; the analysis was run twice to gauge run-to-run variability.

### 6.9 Statistics and software

All *p*-values are two-sided; multiple testing used the Benjamini–Hochberg procedure; rates were compared with Fisher's exact test; correlations are Spearman unless stated. Analyses used R 4.6.1 (`scTenifoldKnk` 1.1, `scTenifoldNet` 1.4, `Matrix` 1.7.5, `locfdr`, `igraph`) and Python 3.10.21 (`scanpy` 1.10.4, `anndata` 0.10.8, `CellOracle` 0.20.0). The stability replicates, the depth-matched sensitivity analysis and the checkpoint evaluation (C3, C4) were run on a second machine with the same R and package versions; only `fgsea` differed (1.38.0 versus 1.34.2), and no enrichment step was re-run.

### 6.10 Data and code availability

All input data are public (GEO: GSE279180; the 10 replication networks are distributed with the `cailab-tamu/scTenifoldKnk` repository). Analysis scripts, the simulation parameters, the 20 stability replicate objects, the complete distance matrices for all knockouts and the CellOracle results accompany this study: code at `[TO BE COMPLETED: repository URL]`, data archived at `[TO BE COMPLETED: Zenodo DOI]`.

### 6.11 Author contributions, funding, competing interests

**Author contributions:** Z.Z. designed the study, performed the analyses and wrote the manuscript. `[TO BE COMPLETED: add co-authors and their contributions if any]`
**Funding:** `[TO BE COMPLETED]`
**Competing interests:** The author declares no competing interests.
**Acknowledgements:** `[TO BE COMPLETED]`

---

## References

All entries below were verified against CrossRef on 2026-09-14 (title, journal, year, volume and DOI).

1. Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Tian Y, Chapkin R, Huang JZ, Cai JJ. scTenifoldKnk: an efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation. *Patterns* 2022;3(3):100434. doi:10.1016/j.patter.2022.100434
2. Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Tian Y, Chapkin R, Huang JZ, Cai JJ. scTenifoldNet: a machine learning workflow for constructing and comparing transcriptome-wide gene regulatory networks from single-cell data. *Patterns* 2020;1(9):100139. doi:10.1016/j.patter.2020.100139
3. Kamimoto K, Stringa B, Hoffmann CM, *et al.* Dissecting cell identity via network inference and in silico gene perturbation. *Nature* 2023;614:742–751. doi:10.1038/s41586-022-05688-9
4. Maggi P, Bulcke CV, Pedrini E, *et al.* B cell depletion therapy does not resolve chronic active multiple sclerosis lesions. *eBioMedicine* 2023;94:104701. doi:10.1016/j.ebiom.2023.104701
5. Wu S, Hu G, Yang Z, Wang Z, Cai J, Mao J, Ge W. Method choice, not biology, determines in silico perturbation results: a systematic evaluation of eight methods across four datasets. Preprint 2026. doi:10.64898/2026.08.11.744106
6. Qiu R, Zhao MM. A confound-diagnostic toolkit for in silico perturbation with single-cell foundation models. Preprint 2026. doi:10.64898/2026.08.04.732812
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B* 1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
8. Efron B. Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. *Journal of the American Statistical Association* 2004;99(465):96–104. doi:10.1198/016214504000000089
9. [Dataset paper for GSE279180] Cell type mapping reveals tissue niches and interactions in subcortical multiple sclerosis lesions. *Nature Neuroscience* 2024. doi:10.1038/s41593-024-01796-z

`[CITATION NEEDED]` — if the reviewers ask for citations on: (i) the original scTenifoldNet application to MS lesions, (ii) alternative in-silico perturbation tools beyond CellOracle, (iii) prior work on weight dilution in network perturbation, please supply or approve a literature search; I did not add unverified references.

---

## Figure legends

**Figure 1. Containment of the differentially regulated set in the knocked-out gene's direct out-neighbours.** (A) |DR| versus outdegree for all 3,195 knockouts in four MS-microglia networks (log–log), with the bound |DR| = outdegree + 1 and the developers' published *Trem2* result (outdegree 2,744; |DR| = 128) as an independent point. (B) Fraction of direct targets called significant versus outdegree, showing weight dilution: ≈1 at outdegree ≤ 10 and 0 for every knockout with outdegree > 10 under the theoretical χ²(1) null. (C) Containment compliance per network (99.8–100.0%; all violations are outdegree-0 knockouts). Unless stated, statistics use the theoretical χ²(1) null and Benjamini–Hochberg control at 5%.

**Figure 2. Replication of containment in 10 published networks and 14 networks in total.** (A) |DR| versus outdegree across the 10 published networks (grey) and the networks of this study (red), with the outdegree + 1 bound. (B) Per-dataset containment rate, 100% in every published network (2,591–9,970 genes; human and mouse; liver, lung, muscle, intestine, brain; knocked-out genes *Hnf4a*+*Hnf4g*, *Nkx2-1*, *Dmd*, *Malat1*, *Ahr*, *Mecp2* ×2, *Cftr* ×2, *Akap7*).

**Figure 3. The containment property is specific to `scTenifoldKnk`.** Comparison with CellOracle on the same microglia data (3,000 cells, 3,011 genes, `hg19_gimmemotifsv5_fpr2`). (A) Number of genes with changed expression after knockout; (B) fraction of affected genes that are direct targets of the knockout (22.1% for *IRF8* to 55.5% for *STAT3*), leaving 44–78% outside the direct target set; (C) run-to-run variability of the two CellOracle runs. Only 4 of the 20 panel genes are transcription factors and therefore eligible for knockout.

**Figure S1. Containment under four parameter settings.** 400-gene Control networks at *nc_lambda* = 0, 0.5, 1 (*nc_q* = 0.90) and *nc_lambda* = 0 (*nc_q* = 0.95); 403 knockouts each. Compliance was 99.0% (4 outdegree-0 violations) under the first three settings and 96.0% (16 outdegree-0 violations) at *q* = 0.95; excluding outdegree-0 knockouts, compliance was 100% in all four settings. Network density and median outdegree change with the settings, the rule does not.

**Figure S2. Stability across independent cell subsamples.** Ten replicates per condition (800 cells sampled per replicate; 10 networks × 500 cells; every gene knocked out). (A) C1 containment per replicate — 100% of non-degenerate knockouts in all 20 replicates (threshold 99%). (B) C2 distribution of the pairwise Spearman ρ of the max-*z* ranking (90 pairs; median 0.185; pre-specified threshold 0.6). (C) C3 smallest candidate-gene FDR per replicate under the empirical null — no candidate gene below 0.05 in any replicate. (D) C4 number of genes at FDR < 0.05 in each replicate-paired Control-versus-MS comparison — 6 of 10 pairs returned 1–2 genes (genes labelled), whereas the primary comparison returned none.

---

## Supplementary tables and data

| Item | Content | File |
|---|---|---|
| Table S1 | Containment results for all 3,195 knockouts | `results/tables/Table_S1_containment_3195_knockouts.csv` |
| Table S2 | Parameter sensitivity (4 configurations, 403 knockouts each) | `results/tables/Table_S2_parameter_sensitivity.csv` |
| Table S3 | Containment in the 10 published networks | `results/tables/Table_S3_replication_10_datasets.csv` |
| Tables S4, S5 | Empirical *p*-values for the panel (Control, MS) | `results/tables/Table_S4_…`, `Table_S5_…` |
| Table S6 | Specificity metrics (max-*z*, *n*~z>2~) per gene and condition | `results/tables/Table_S6_specificity_metric.csv` |
| Tables S7, S8 | Panel summary under the empirical null (Control, MS) | `results/tables/Table_S7_…`, `Table_S8_…` |
| Table S9 | CellOracle affected genes versus direct targets | `results/tables/Table_S9_celloracle_vs_direct_targets.csv` |
| Table S10 | Per-replicate stability metrics (20 replicates) | `results/tables/Table_S10_stability_replicates.csv` |
| Table S11 | Stability checkpoint verdicts C1–C4 | `results/tables/Table_S11_stability_summary.csv` |
| Table S12 | C3 checkpoint: empirical-null panel result per replicate | `outputs/results/vko_stability_checkpoints/checkpoint_C3_empirical_null_by_replicate.csv` |
| Table S13 | C4 checkpoint: paired condition comparison per replicate | `outputs/results/vko_stability_checkpoints/checkpoint_C4_paired_condition_comparison.csv` |
| Supplementary Data 1 | Depth-matched versus primary condition comparison | `outputs/results/vko_ms_depthmatched/panel_depthmatched_vs_main.csv` |

Figures: `results/figures/Fig_degree_artifact.{png,pdf}` (Figure 1), `Fig_replication.*` (Figure 2), `Fig_tool_comparison.*` (Figure 3), `FigS1_parameter_sensitivity.*`, `FigS2_stability.*`.

---

## Internal appendix (not for submission): claim → evidence map

| Statement in the manuscript | Source file |
|---|---|
| 3,191 / 3,195 knockouts contained (99.87%); 4 violations, all outdegree 0 | `Table_S1_containment_3195_knockouts.csv` |
| 3,064 knockouts returned exactly 1 significant gene under the χ² null; max 15 | same |
| `frac_targets_sig` = 0.92 for outdegree 1–10 and 0 for outdegree > 10 | same (computed) |
| Panel under empirical null: `BTK` 197→198, `CD33` 146→147, `YKT6` 145→146; hubs 767→1; `ACTA2` (MS, outdegree 0) → 188 genes (adj. *p* 1.3 × 10⁻²¹⁶) | `Table_S7`, `Table_S8`, `MSvC_MS_outdegree.csv` |
| χ²-null outdegree-0 knockouts returned 7–15 genes (TSHR 15, ACTA2 7, FAF1 7, TSHR 9) | `Table_S1_…csv` |
| Published *Trem2* result: outdegree 2,744, |DR| = 128, 128/128 inside direct targets | `docs/06_core_finding_CN.md`, `code/04_replication/33_author_trem2_validation.R` |
| Replication: 100% containment in 10 published networks; `n_DR` 2–377 | `Table_S3_replication_10_datasets.csv` |
| Parameter sensitivity: 99.0% / 96.0% compliance; 100% after excluding outdegree-0 | `Table_S2_parameter_sensitivity.csv` |
| Condition comparison 7.25% versus split-half null 5.65%; Fisher *p* = 0.221, OR 1.31 (0.86–2.00); Jaccard 0.136 | `docs/05_condition_comparison_CN.md` |
| No panel gene at FDR < 0.05 under the empirical null (min adj. *p* 0.525, Control) | `Table_S4`, `Table_S5` |
| CellOracle: 22.1% (*IRF8*) to 55.5% (*STAT3*) of affected genes are direct targets | `Table_S9_celloracle_vs_direct_targets.csv` |
| Stability: C1 100% (20/20), C2 median ρ 0.185 (90 pairs), C3 0/20 (min candidate FDR 0.062), C4 6/10 pairs | `Table_S10`, `Table_S11`, checkpoint CSVs |
| Depth matching: ratio 2.07 → ≈1.20 (post-QC; 3.84 → 1.20 over all cells); 0/800 at FDR < 0.05; network ρ 0.317 versus 0.196; panel median *p* 0.454 → 0.595 | `outputs/results/vko_ms_depthmatched/panel_depthmatched_vs_main.csv`, `docs/13_…`, `vko-run/depth_window_numbers.R` |
| Data: 23,627 genes × 9,239 microglia; 2,005 Control / 7,234 MS; median 927 / 3,564 counts; 1,347 / 6,927 cells after QC | `docs/02_data_landing_CN.md`, `work/data/metadata.csv`, `vko-run/verify_qc_numbers.R` |
