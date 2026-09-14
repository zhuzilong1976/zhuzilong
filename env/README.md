# Recorded software environment

| File | Contents |
|---|---|
| `R_sessionInfo.txt` | Full `sessionInfo()` of the R 4.6.1 installation used for all statistical analyses |
| `python_requirements.txt` | `pip freeze` of the isolated CPython 3.10.21 environment used for the CellOracle comparison (174 packages) |
| `node_version.txt` | Node.js version (v24.19.0) used for the download helpers |

## Key versions

| Software | Version |
|---|---|
| R | 4.6.1 (2026-06-24 ucrt), x86_64-w64-mingw32 |
| scTenifoldKnk | 1.1 (CRAN, published 2026-09-02) |
| scTenifoldNet | 1.4 |
| Matrix | 1.7.5 |
| locfdr | 1.1-8 |
| enrichR | 3.4 |
| igraph | 2.3.3 |
| fgsea | 1.34.2 (Bioconductor) |
| Python | 3.10.21 |
| scanpy | 1.10.4 |
| anndata | 0.10.8 |
| numpy | 1.26.4 |
| pandas | 1.5.3 |
| CellOracle | 0.20.0 |

## Platform notes

All analyses were run on Windows 10 (build 19045), 8 cores, without
administrator rights.

Two environment constraints shaped the workflow and should be noted when
reproducing:

1. **R could not use TLS in this sandbox** (schannel credential error), so CRAN
   binary packages were downloaded with Node.js and installed offline
   (`code/00_setup/fetch_r_packages.mjs`, `install_r_packages_offline.R`).
   On a normal machine, `install.packages("scTenifoldKnk")` plus the packages
   listed above is sufficient.
2. **CellOracle has no Windows-compatible dependency set** (`velocyto`,
   `gimmemotifs`, `pybedtools` are source-only and require a C compiler), so it
   was installed into a portable CPython 3.10 with three placeholder modules.
   The placeholders raise on use and were never invoked; see
   `docs/07_tool_comparison_CN.md` and limitation 7 of
   `docs/09_methods_and_limitations_draft_EN.md`.
