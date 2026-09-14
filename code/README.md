# Code

Scripts are numbered in suggested execution order. Unless stated otherwise they
are run **from the repository root** (`Rscript code/<dir>/<script>.R`) and write
their outputs under `results/` or `work/`.

## Path conventions

The scripts are published **as they were executed**, and therefore use the
original directory convention of the analysis project:

| Directory | Role |
|---|---|
| `work/` | inputs, downloaded resources and intermediate objects (created on first run) |
| `outputs/results/` | analysis outputs written by the scripts |

In this repository the curated final products are additionally exposed at the
root as `results/figures/` and `results/tables/` (a copy of the corresponding
files from `outputs/results/`). No script paths were rewritten, so that the code
that produced the published numbers is the code that is distributed; when
re-running, create a working directory that contains `work/` and `outputs/` and
execute the scripts from there, or adjust the two constants below.

The only paths that need changing for a fresh checkout are the input locations
declared at the top of each script, which follow the pattern:

```r
counts_file <- "work/data/counts.mtx"          # produced by 01_data/01_…
out_dir     <- "outputs/results/vko_..."       # where results are written
```

and, for the Python scripts, `work/celloracle/…` for the prepared AnnData object.

Column "Run" gives the observed runtime on an 8-core Windows 10 machine; network
construction is effectively single-threaded in the package as used here.

---

## 00_setup — environment

| Script | Purpose | Run |
|---|---|---|
| `00_setup_environment.R` | Checks R/package versions and runs a small end-to-end smoke test (simulated counts, 100 genes × 500 cells). Prints a timing baseline and writes `results/tables/environment_report.txt`. | < 1 min |
| `fetch_r_packages.mjs` | Resolves the CRAN dependency closure for a set of packages and downloads Windows binary `.zip` files. Used because the analysis machine had no working TLS stack for R. | minutes |
| `install_r_packages_offline.R` | Installs the downloaded `.zip` files in topological order. | minutes |
| `resolve_python_wheels.mjs` | Resolves and downloads a Python dependency closure compatible with CPython 3.10 / win_amd64 (used to build the CellOracle environment). | minutes |

For a standard installation with internet access, the first script alone is
sufficient; the three offline helpers are only needed to reproduce the exact
environment used here.

## 01_data — input data and cohort diagnostics

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `01_prepare_microglia_data.R` | Converts the GEO microglia AnnData into a counts matrix plus metadata. | `data/GSE279180_ctype_MG.h5ad` | `work/data/{counts.mtx,genes.txt,barcodes.txt,metadata.csv}` |
| `02_check_panel_expression.R` | Expression rate of every panel gene overall and by condition; identifies genes that would be removed by the QC filter. | `work/data/*` | console |
| `03_select_negative_controls.R` | Evaluates candidate negative-control genes against the expression filter. | `work/data/*` | console |
| `04_negative_control_search.R` | Data-driven search for low-variance, moderately expressed genes. | `work/data/*` | console |
| `05_library_size_by_condition.R` | Library-size distributions per condition and per sample; retention under alternative thresholds. | `work/data/*` | console |
| `06_check_batch_confounding.R` | Tests whether sequencing batch is confounded with condition. | `work/data/*` | console |
| `07_check_donor_confounding.R` | Tests whether lesion type / condition is confounded with donor. **This analysis is why the chronic-active vs chronic-inactive design was rejected.** | `work/data/*` | console |

## 02_network_analysis — network construction, virtual knockout, containment

| Script | Purpose | Run |
|---|---|---|
| `10_run_vko_panel.R` | **Main pipeline.** QC → shared gene set → network construction → virtual knockout of a gene panel or of a chosen subset, in transcriptome-wide mode (one network, many knockouts). Writes distance matrices, per-gene statistics, outdegree, negative-control diagnostics and a run log. Accepts a configuration file as its first argument (see `code/config/`). | ≈2 h per condition (800 genes, 10 networks × 500 cells) |
| `11_containment_all_genes.R` | **Core test.** For every gene in a network, knocks it out, recomputes the differential-regulation statistics from the saved distance matrix, and tests DR ⊆ {KO} ∪ direct targets and \|DR\| ≤ outdegree + 1. | ≈10 s per network |
| `12_containment_panel_genes.R` | The same test for the 20 panel genes, with a per-gene printout (used for the first observation of the rule). | seconds |
| `13_outdegree_diagnostic.R` | Reports outdegree for panel genes and flags genes with zero outdegree, whose knockout is a no-op. | seconds |
| `14_degree_artifact.R` | Quantifies the relationship between the number of significant genes and outdegree (the "= outdegree + 1" pattern). | seconds |
| `15_diagnose_propagation.R` | Checks whether the knockout signal propagates beyond the knocked-out gene (max/median FC per knockout). | seconds |
| `16_test_propagation_vs_density.R` | Tests whether network sparsification restores propagation (it does not). | minutes |

## 03_controls — null models and sensitivity analyses

| Script | Purpose | Run |
|---|---|---|
| `20_full_knockout_distances.R` | Knocks out **every** gene of a network and stores the full *n* × *n* distance matrix. Usage: `Rscript … <tag>`; takes a network `.rds` and an output path. Later scripts reuse these matrices, so this only runs once per network. | ≈7 min per network (800 genes) |
| `21_empirical_null_analysis.R` | Builds the empirical null for max_z from all non-panel knockouts and computes empirical *p*-values (overall and outdegree-matched, ±20%) for the panel genes. | seconds |
| `22_build_split_half_networks.R` | Builds a network from a random half of the Control cells. Usage: `<tag> <condition> <A\|B>`. | ≈2 h per network |
| `23_split_half_null_comparison.R` | **Technical null.** Compares two same-condition half-sample networks and contrasts the significant-gene rate with the disease-vs-control comparison. | seconds |
| `24_condition_comparison.R` | Two-condition comparison (Control network vs MS network) via manifold alignment and differential regulation — the framework's original design use. | seconds (networks already built) |
| `25_parameter_sensitivity.R` | Rebuilds 400-gene networks under four parameter settings and repeats the containment test. Usage: `<tag> <lambda> <q>`. | ≈21 min per configuration |
| `26_scan_depth_windows.R` | Scans candidate library-size windows for the trade-off between cell yield and depth matching. | seconds |
| `27_depth_matched_sensitivity.R` | Repeats the condition comparison using depth-matched networks (window 1,500–5,000 counts). | seconds (after network builds) |

## 04_replication — cross-dataset replication

| Script | Purpose |
|---|---|
| `30_fetch_published_networks.mjs` | Downloads the 10 network objects published with scTenifoldKnk into `work/author_networks/` (108 MB). |
| `31_inspect_published_networks.R` | Lists the objects and structure of each downloaded `.RData` file. |
| `32_replication_containment.R` | **Core replication test.** Identifies the knocked-out gene from the data (rows non-zero in `WT` but all-zero in `KO`) and tests containment for each dataset. |
| `33_author_trem2_validation.R` | Inspects the developers' published TREM2 result object (network, manifold alignment, DR table). |
| `34_fetch_author_results.mjs` | Downloads the published differential-regulation tables (used to compare significant-gene counts across tools). |

## 05_tool_comparison — CellOracle (Python)

Requires the isolated CPython 3.10 environment described in
`docs/07_tool_comparison_CN.md`. All three scripts must be run with that
interpreter.

| Script | Purpose | Run |
|---|---|---|
| `40_celloracle_prepare_data.py` | Prepares the AnnData object for CellOracle: intersect with the base GRN genes, normalise, log-transform, select 3,000 HVGs plus all panel genes. | ≈1 min |
| `41_celloracle_run.py` | Builds the Oracle object, imputes expression, constructs the GRN (`bagging_number = 5`), fits it for simulation, and simulates knockout of every panel gene that is a transcription factor (`n_propagation = 3`). Saves the GRN edge table and the per-cell expression-change matrices. | ≈12 min |
| `42_celloracle_compare.py` | Compares the affected-gene sets with the TF's direct targets in the inferred GRN. | seconds |

## 06_figures

| Script | Figure | Notes |
|---|---|---|
| `50_figure_degree_containment.R` | Figure 1 | Recomputes containment for the four MS-microglia networks, then draws the three panels; also writes `results/figures/containment_all_networks.csv` (3,195 knockouts). |
| `51_figure_parameter_sensitivity.R` | Figure S1 | Reads the four `code/03_controls/25_…` outputs. |
| `52_figure_tool_comparison.R` | Figure 3 | Reads `Table_S9`; the run-to-run values for the second CellOracle run are hard-coded in the script (both runs documented in `docs/07`). |
| `53_figure_cross_dataset.R` | Figure 2 | Reads `Table_S1` and `Table_S3`. |
| `54_figure_stability.R` | Figure S2 | Four-panel summary of the stability checkpoints C1–C4 from `work/stability/*_rep*.rds`, `Table_S10` and (for C3/C4) the checkpoint files in `outputs/results/vko_stability_checkpoints/`, which are produced by `07_supplementary/67_stability_checkpoints_C3C4.R`. |

All figures are written as 2,400-px PNG and vector PDF, with English labels.

## 07_supplementary

| Script | Purpose |
|---|---|
| `60_perturbation_profile_comparison.R` | Correlates perturbation profiles of the same knockout between the Control and MS networks. |
| `61_target_enrichment.R` | Hypergeometric enrichment of each knockout's top-ranked targets against KEGG/Reactome/GO/Hallmark (background = network genes). |
| `62_target_enrichment_specificity.R` | The same after per-target specificity normalisation (z across knockouts), which is required because the raw ranking is dominated by a common component. |
| `63_specificity_metric.R` | Computes max_z and n_{z>2} per knockout. |
| `64_top_targets_per_knockout.R` | Lists the most specific targets per knockout. |
| `65_fisher_test_null_vs_observed.R` | Fisher's exact test comparing the split-half null with the observed comparison. |
| `66_summary_tables.R` | Assembles the cross-condition summary table. |
| `67_stability_checkpoints_C3C4.R` | Completes the two stability checkpoints that `03_controls/29_stability_analysis.R` leaves as pointers: C3 (panel genes under the empirical null, per replicate) and C4 (paired Control-vs-MS condition comparison per replicate). Reads `work/stability/*_rep*.rds` and writes the two checkpoint CSVs consumed by `06_figures/54_figure_stability.R`. |

## config

`10_run_vko_panel.R` sources a configuration file that overrides the defaults:

| File | Network | Cells | Window |
|---|---|---|---|
| `config_control.R` | MS-microglia Control | all passing QC (minLib 500) | — |
| `config_MS.R` | MS-microglia MS | all passing QC (minLib 500) | — |
| `config_control_depthmatched.R` | MS-microglia Control | depth-matched | 1,500–5,000 counts |
| `config_MS_depthmatched.R` | MS-microglia MS | depth-matched | 1,500–5,000 counts |
| `config_smoke_test.R` | reduced smoke test | 200 cells | — |

## utils

Small helpers used during development (Enrichr gene-set download, GEO
supplementary-file listing, CSV validation) are kept in `code/utils/` for
completeness. They are not required to reproduce the reported results.
