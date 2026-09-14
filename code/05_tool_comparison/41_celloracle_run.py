"""Run CellOracle GRN construction + in-silico KO on the MS microglia data."""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
from pathlib import Path
import sys

OUT = Path("outputs/results/celloracle")
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ["TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
         "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA"]

print("=== 1. 载入预处理数据 ===", flush=True)
adata = sc.read_h5ad("work/celloracle/GSE279180_MG_hvg.h5ad")
print("数据:", adata.shape, flush=True)
# 单一小胶质群体作为一个 GRN 单元（CellOracle 教程的标准用法；
# 亚型特异 GRN 需要显著更多计算时间，且不影响我们这里的比较问题）
adata.obs["cell_type"] = "Microglia"
# 计算预算限制：下采样到 3,000 细胞（CellOracle 教程常用规模）
N_CELLS = 3000
if adata.n_obs > N_CELLS:
    import numpy as _np
    _rng = _np.random.default_rng(1)
    # 保持 MS / Control 比例
    idx = []
    for cond, grp in adata.obs.groupby("condition"):
        n = int(round(N_CELLS * len(grp) / adata.n_obs))
        idx.extend(_rng.choice(grp.index.values, size=min(n, len(grp)), replace=False))
    adata = adata[_np.array(idx)].copy()
    print(f"下采样后: {adata.shape}", flush=True)

print("\n=== 2. 构建 Oracle 对象 ===", flush=True)
base_grn = co.data.load_human_promoter_base_GRN()
tf_names = set(base_grn.columns) - {"peak_id", "gene_short_name"}
print(f"base GRN 中的转录因子数: {len(tf_names)}", flush=True)
in_tf = [g for g in PANEL if g in tf_names]
not_tf = [g for g in PANEL if g not in tf_names]
print(f"面板中属于 TF（CellOracle 可敲除）: {in_tf}", flush=True)
print(f"面板中不属于 TF（CellOracle 无法敲除）: {not_tf}", flush=True)

oracle = co.Oracle()
oracle.import_anndata_as_normalized_count(
    adata=adata, cluster_column_name="cell_type", embedding_name="X_umap"
)
print("表达数据已导入；细胞数:", oracle.adata.shape, flush=True)
oracle.import_TF_data(TF_info_matrix=base_grn)
print("TF 数据已导入", flush=True)

print("\n=== 3. 插补与 PCA ===", flush=True)
oracle.perform_PCA()
oracle.knn_imputation(n_pca_dims=50, k=30, n_jobs=-1)
print("插补完成", flush=True)

print("\n=== 4. 构建 GRN ===", flush=True)
links = oracle.get_links(cluster_name_for_GRN_unit="cell_type", alpha=10,
                         bagging_number=5, verbose_level=1, n_jobs=-1)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=10000)
print("GRN 完成", flush=True)

# 保存 GRN（TF -> target 的边表），用于计算直接靶标
try:
    ld = links.links_dict
    for cl, df in ld.items():
        df.to_csv(OUT / f"GRN_links_{cl}.csv", index=False)
        print(f"  GRN 边表已保存: {cl} -> {df.shape}", flush=True)
except Exception as e:
    print("  GRN 边表保存失败:", e, flush=True)
try:
    oracle.to_hdf5(OUT / "oracle_object.hdf5")
    print("  Oracle 对象已保存", flush=True)
except Exception as e:
    print("  Oracle 保存失败（不影响结果）:", e, flush=True)

print("\n=== 5. 为模拟准备 GRN ===", flush=True)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
print("模拟用 GRN 就绪", flush=True)

print("\n=== 6. 逐个敲除 ===", flush=True)
results = {}
for g in in_tf:
    try:
        oracle.simulate_shift(perturb_condition={g: 0.0}, n_propagation=3)
        avail = list(oracle.adata.layers.keys())
        if "delta_X" not in avail:
            raise KeyError(f"delta_X 缺失；可用层: {avail}")
        delta = np.asarray(oracle.adata.layers["delta_X"]).copy()
        results[g] = {
            "delta": delta,
            "genes": list(oracle.adata.var_names),
        }
        n_aff = int((np.abs(delta).sum(axis=0) > 0).sum())
        print(f"  {g}: 完成；delta_X 非零基因数 {n_aff}；层: {avail}", flush=True)
    except Exception as e:
        print(f"  {g}: 失败 {type(e).__name__}: {e}", flush=True)

import pickle
with open(OUT / "celloracle_ko_results.pkl", "wb") as f:
    pickle.dump({"results": results, "in_tf": in_tf, "not_tf": not_tf,
                 "var_names": list(oracle.adata.var_names)}, f)
print("\n已保存:", OUT / "celloracle_ko_results.pkl", flush=True)
