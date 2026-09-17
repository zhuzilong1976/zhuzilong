# Containment of virtual-knockout predictions in single-cell gene regulatory networks

Code and derived data accompanying a methodological study of *in-silico* gene
knockout on single-cell gene regulatory networks (GRNs), using human multiple
sclerosis (MS) lesion microglia as the worked example.

## Summary of the finding

The "differentially regulated" (DR) gene set produced by single-cell virtual
knockout is **not a biological readout of downstream consequences**: it is
contained in the knocked-out gene's **direct out-neighbours** in the inferred
network.

> **DR(*g*) ⊆ {*g*} ∪ {direct out-neighbours of *g*}**, equivalently |DR(*g*)| ≤ outdegree(*g*) + 1

This containment held in **38 networks and ~20,800 knockouts** (99.9%), spanning
human and mouse, five tissues, network sizes of 800–9,970 genes, and outdegrees
from 1 to 3,632. This count includes 20 cell-subsample replicates built to test
stability (10 per condition). The only exceptions were genes with zero outdegree,
where the knockout is a no-op and the test statistic degenerates.

Whether a target is called significant is governed by **weight dilution**: a
target passes the significance threshold in proportion to the fraction of its
total incoming weight that the knocked-out gene contributes (Spearman-consistent
monotone decrease of |DR|/outdegree with outdegree, from ≈1 at outdegree 1 to
0.05 at outdegree 3,632).

The behaviour is **tool-specific, not generic to GRN-based knockout simulation**:
with CellOracle (a different design that propagates the perturbation
iteratively), only 22–56% of affected genes are direct targets of the knocked-out
TF.

Three further results matter for practice:

1. **A same-condition split-half comparison produces as many "significant" genes
   as a disease-vs-control comparison** (5.65% vs 7.25%; Fisher *p* = 0.22), and
   the top genes of the null comparison are as biologically plausible as those of
   the disease comparison. Plausibility of a gene list is therefore not evidence.
2. **The candidate-gene panel was indistinguishable from random network genes**
   under an empirical null built from all knockouts (no gene significant after
   FDR correction in either network).
3. **Stability is property-dependent.** Across 10 independent cell subsamples per
   condition, the qualitative results replicated exactly — containment held for
   100% of non-degenerate knockouts in all 20 replicates, and no panel gene was
   significant under the empirical null in any replicate. Magnitude-based outputs
   did not replicate: the median correlation of the per-gene perturbation ranking
   between replicates was 0.185, and 6 of 10 replicate-paired condition
   comparisons produced 1–2 genes at FDR < 0.05 while the primary comparison
   produced none (Table S10/S11, Figure S2).

## Repository layout

```
.
├── README.md                     this file
├── LICENSE                       MIT (code); CC-BY-4.0 for figures/tables
├── CITATION.cff
├── code/                         all analysis code (see code/README.md)
│   ├── 00_setup/                 environment construction and offline installers
│   ├── 01_data/                  data preparation and cohort diagnostics
│   ├── 02_network_analysis/      network construction, virtual knockout, containment
│   ├── 03_controls/              split-half null, empirical null, condition comparison,
│   │                             parameter and depth-matched sensitivity
│   ├── 04_replication/           cross-dataset replication on published networks
│   ├── 05_tool_comparison/       CellOracle comparison (Python)
│   ├── 06_figures/               figure generation
│   ├── 07_supplementary/         enrichment, specificity metrics, summary tables
│   ├── config/                   analysis parameter files
│   └── utils/                    download helpers
├── docs/                         analysis notes and the Methods/Limitations draft
├── env/                          recorded software versions
├── data/                         input data (not tracked; see data/README.md)
└── results/
    ├── figures/                  Figures 1–3 and Supplementary Figures S1–S2 (PNG + vector PDF)
    └── tables/                   Supplementary Tables S1–S11
```

## Data

All input data are public. See `data/README.md` for download instructions.
In brief:

* **Primary dataset** — `GSE279180`, human subcortical MS lesions, microglia
  subset (`GSE279180_ctype_MG.h5ad`, 23,627 genes × 9,239 cells, raw counts).
* **Replication networks** — 10 network objects published by the developers of
  scTenifoldKnk in `inst/manuscript/*/Results/*.RData` of
  `cailab-tamu/scTenifoldKnk`.
* **CellOracle base GRN** — `hg19_gimmemotifsv5_fpr2` (downloaded automatically
  by the `celloracle` package).

## Reproducing the analysis

The pipeline is designed to be run from the repository root. Scripts are
numbered in execution order; each one documents its inputs and outputs in its
header and is listed in `code/README.md`.

```bash
# 1. environment
Rscript code/00_setup/00_setup_environment.R

# 2. data preparation (requires data/GSE279180_ctype_MG.h5ad)
Rscript code/01_data/01_prepare_microglia_data.R
Rscript code/01_data/02_check_panel_expression.R
Rscript code/01_data/03_select_negative_controls.R

# 3. primary network analysis (≈2 h per condition on 8 cores)
Rscript code/02_network_analysis/10_run_vko_panel.R code/config/config_control.R
Rscript code/02_network_analysis/10_run_vko_panel.R code/config/config_MS.R

# 4. containment and diagnostics
Rscript code/02_network_analysis/11_containment_all_genes.R
Rscript code/02_network_analysis/14_degree_artifact.R

# 5. controls
Rscript code/03_controls/20_full_knockout_distances.R MSvC_Control
Rscript code/03_controls/21_empirical_null_analysis.R
Rscript code/03_controls/22_build_split_half_networks.R CtrlHalfA Control A
Rscript code/03_controls/23_split_half_null_comparison.R
Rscript code/03_controls/24_condition_comparison.R

# 5b. depth-matched sensitivity (window 1,500-5,000 counts, ~2 h per condition)
Rscript code/02_network_analysis/10_run_vko_panel.R code/config/config_control_depthmatched.R
Rscript code/02_network_analysis/10_run_vko_panel.R code/config/config_MS_depthmatched.R
Rscript code/03_controls/27_depth_matched_sensitivity.R

# 5c. stability: 10 independent cell subsamples per condition (~2-4 h each;
#     launchers: run_stability_local.ps1 / run_stability_linux.sh / run_stability_slurm.sh)
Rscript code/03_controls/28_stability_replicates.R MSvC_Control 1 800 800
Rscript code/03_controls/29_stability_analysis.R              # Table S10, S11 (C1, C2)
Rscript code/07_supplementary/67_stability_checkpoints_C3C4.R # checkpoints C3, C4

# 6. replication on published networks
node code/04_replication/30_fetch_published_networks.mjs
Rscript code/04_replication/32_replication_containment.R

# 7. tool comparison (Python; see code/README.md for the environment)
python code/05_tool_comparison/40_celloracle_prepare_data.py
python code/05_tool_comparison/41_celloracle_run.py
python code/05_tool_comparison/42_celloracle_compare.py

# 8. figures
Rscript code/06_figures/50_figure_degree_containment.R
Rscript code/06_figures/51_figure_parameter_sensitivity.R
Rscript code/06_figures/52_figure_tool_comparison.R
Rscript code/06_figures/53_figure_cross_dataset.R
Rscript code/06_figures/54_figure_stability.R
```

Scripts write their outputs to `outputs/results/` and their intermediates to
`work/`, which is created on first run; both are relative to the directory the
scripts are executed from. The curated final products are also exposed at the
repository root as `results/figures/` and `results/tables/`. See the "Path
conventions" section of `code/README.md` for details — the scripts are
distributed exactly as they were executed, with no paths rewritten.

## Software

Recorded versions are in `env/`. Principal packages:

| Software | Version | Purpose |
|---|---|---|
| R | 4.6.1 | all statistical analysis |
| scTenifoldKnk | 1.1 (CRAN, 2026-09-02) | network construction, virtual knockout, differential regulation |
| scTenifoldNet | 1.4 | network construction, manifold alignment |
| locfdr | 1.1-8 | Efron empirical null |
| fgsea | 1.34.2 | pathway enrichment |
| Python | 3.10.21 | CellOracle comparison |
| CellOracle | 0.20.0 | second-method comparison |

The R environment used here was installed without administrator rights by
extracting the official Windows installer and installing 62 CRAN binary
packages offline (`code/00_setup/`). The CellOracle environment used a portable
CPython distribution plus three placeholder modules for packages with no Windows
build (`velocyto`, `gimmemotifs`, `pybedtools`); the placeholders raise on use and
were never invoked by the code paths reported here — see
`docs/07_tool_comparison_CN.md` for details and `docs/09_methods_and_limitations_draft_EN.md`
for the associated limitation.

## Known gaps

* The primary analysis uses one dataset; the candidate panel was not enriched in
  the primary networks.
* The depth-matched sensitivity analysis (window 1,500–5,000 counts) lowers the
  Control/MS median-depth ratio (2.07 measured after quality control; 3.84 → 1.20
  over all cells) but does not fully equalise the two depth distributions.
* The stability replicates are smaller than the primary networks (343–369 vs
  1,347 Control cells after quality control), so part of the variability they
  measure reflects sample size rather than a property of the primary networks.

## Citation

See `CITATION.cff`. If you use the containment analysis or the figures, please
cite this repository together with the scTenifoldKnk paper (Osorio *et al.*,
*Patterns* 2022;3:100434).

Repository: <https://github.com/zhuzilong1976/zhuzilong>
Archived release: <https://doi.org/10.5281/zenodo.22752977> (concept DOI, always
the latest version; v1.1.0 is <https://doi.org/10.5281/zenodo.22814877> and v1.0.0
is <https://doi.org/10.5281/zenodo.22752978>)
