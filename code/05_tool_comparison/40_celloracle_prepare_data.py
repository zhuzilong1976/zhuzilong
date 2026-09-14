"""Prepare MS-microglia AnnData for CellOracle using the recommended gene budget.

CellOracle is designed for ~1,000-3,000 genes. We therefore keep the highly
variable genes (plus all panel genes present in the base GRN), which is the
standard usage shown in the CellOracle tutorials.
"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
from pathlib import Path

OUT = Path("work/celloracle")
OUT.mkdir(parents=True, exist_ok=True)
N_HVG = 3000
PANEL = ["TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
         "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA"]

print("=== 1. 读入原始数据 ===", flush=True)
adata = sc.read_h5ad("work/geo/GSE279180_ctype_MG.h5ad")
adata.X = adata.X.tocsr()
adata.layers["raw_count"] = adata.X.copy()

print("=== 2. 与 base GRN 基因取交集 ===", flush=True)
base_grn = co.data.load_human_promoter_base_GRN()
grn_genes = set(base_grn["gene_short_name"].astype(str))
adata = adata[:, adata.var_names.isin(grn_genes)].copy()
print("交集后:", adata.shape, flush=True)

print("=== 3. 预处理 ===", flush=True)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
hvg = adata.var["highly_variable"].values
keep = hvg | adata.var_names.isin(PANEL)
print(f"高变基因 {int(hvg.sum())}，并上面板基因后保留 {int(keep.sum())}", flush=True)
adata = adata[:, keep].copy()

sc.tl.pca(adata, n_comps=50, use_highly_variable=False)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
adata.obs["cell_type"] = adata.obs["subtype"].astype(str)

panel_in = [g for g in PANEL if g in adata.var_names]
missing = [g for g in PANEL if g not in adata.var_names]
print(f"面板基因保留 {len(panel_in)}/20；缺失: {missing}", flush=True)
print("最终数据:", adata.shape, flush=True)

adata.write_h5ad(OUT / "GSE279180_MG_hvg.h5ad")
print("已保存:", OUT / "GSE279180_MG_hvg.h5ad", flush=True)
