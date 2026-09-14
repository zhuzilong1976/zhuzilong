# GSE279180 小胶质细胞 h5ad -> scTenifoldKnk 输入格式
# 输出：work/data/{counts.mtx, genes.txt, barcodes.txt, metadata.csv}

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
})

h5_path <- "work/geo/GSE279180_ctype_MG.h5ad"
out_dir <- "work/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

f <- H5File$new(h5_path, mode = "r")

## ---- 1. 基因名与细胞名 ----
genes <- f[["var/_index"]]$read()
barcodes <- f[["obs/_index"]]$read()
n_genes <- length(genes)
n_cells <- length(barcodes)
cat(sprintf("基因数 = %d，细胞数 = %d\n", n_genes, n_cells))
cat("基因名前 8 个：", paste(head(genes, 8), collapse = ", "), "\n")

## ---- 2. 读稀疏矩阵（AnnData: X 为 cells x genes 的 CSR）----
ptr <- f[["X/indptr"]]$read()
idx <- f[["X/indices"]]$read()
val <- f[["X/data"]]$read()
stopifnot(length(ptr) == n_cells + 1)
cat(sprintf("非零元素 = %s (%.1f%% 稀疏度)\n",
            format(length(val), big.mark = ","),
            100 * (1 - length(val) / (n_cells * n_genes))))

rows <- rep(seq_len(n_cells), times = diff(ptr))
X_cells_by_genes <- sparseMatrix(
  i = rows, j = idx + 1L, x = val,
  dims = c(n_cells, n_genes),
  dimnames = list(barcodes, genes)
)
rm(rows, idx, val, ptr)

## scTenifoldKnk 需要：基因(行) x 细胞(列)
X <- t(X_cells_by_genes)
X <- as(X, "CsparseMatrix")
rm(X_cells_by_genes)
gc()

cat(sprintf("转换为 %d 基因 x %d 细胞，非零 %s\n", nrow(X), ncol(X),
            format(length(X@x), big.mark = ",")))
cat(sprintf("数值范围：%g – %g；是否全为整数：%s\n",
            min(X@x), max(X@x), all(X@x == round(X@x))))

## ---- 3. 元数据 ----
read_obs <- function(f, name) {
  o <- f[["obs"]][[name]]
  if (inherits(o, "H5Group")) {
    cats <- o[["categories"]]$read()
    codes <- as.integer(o[["codes"]]$read())
    out <- ifelse(codes < 0, NA, cats[codes + 1L])
    return(out)
  }
  o$read()
}

meta <- data.frame(
  barcode     = barcodes,
  celltype    = "Microglia",
  condition   = read_obs(f, "condition"),
  donor       = read_obs(f, "patient_id"),          # 样本级（患者+样本编号）
  patient     = NA_character_,                      # 患者级（去掉样本编号）
  lesion_type = read_obs(f, "lesion_type"),
  subtype     = read_obs(f, "subtype"),
  sample_id   = read_obs(f, "sample_id"),
  age         = read_obs(f, "age"),
  sex         = read_obs(f, "sex"),
  batch_sn    = read_obs(f, "batch_sn"),
  stringsAsFactors = FALSE
)
## 患者级 ID：patient_id 形如 "MS377 A2D4"，去掉空格后的样本编号
meta$patient <- sub("[[:space:]].*$", "", meta$donor)
f$close_all()

## ---- 4. 写出 ----
Matrix::writeMM(X, file.path(out_dir, "counts.mtx"))
writeLines(rownames(X), file.path(out_dir, "genes.txt"))
writeLines(colnames(X), file.path(out_dir, "barcodes.txt"))
utils::write.csv(meta, file.path(out_dir, "metadata.csv"), row.names = FALSE)

## ---- 5. 汇总 ----
cat("\n===== 元数据分布 =====\n")
cat("\n[condition]\n"); print(table(meta$condition, useNA = "ifany"))
cat("\n[lesion_type]\n"); print(table(meta$lesion_type, useNA = "ifany"))
cat("\n[subtype]\n"); print(table(meta$subtype, useNA = "ifany"))
cat("\n[condition x lesion_type]\n"); print(table(meta$condition, meta$lesion_type))
cat("\n[每个供体的细胞数]\n"); print(sort(table(meta$donor), decreasing = TRUE))
cat("\n[每个患者的细胞数]\n"); print(sort(table(meta$patient), decreasing = TRUE))
cat("\n[MS 组内患者占比]\n")
ms_pat <- sort(table(meta$patient[meta$condition == "MS"]), decreasing = TRUE)
print(round(100 * ms_pat / sum(ms_pat), 1))

cat("\n===== 关键门槛检查 =====\n")
ctrl <- sum(meta$condition == "Control")
ms <- sum(meta$condition == "MS")
cat(sprintf("Control 小胶质: %d 个细胞 (>=500: %s)\n", ctrl, ctrl >= 500))
cat(sprintf("MS 小胶质:      %d 个细胞 (>=500: %s)\n", ms, ms >= 500))
cat(sprintf("CA(慢性活动病灶): %d 个细胞\n", sum(meta$lesion_type == "CA", na.rm = TRUE)))
cat(sprintf("CI(慢性非活动):   %d 个细胞\n", sum(meta$lesion_type == "CI", na.rm = TRUE)))
cat(sprintf("Ctrl(对照组织):   %d 个细胞\n", sum(meta$lesion_type == "Ctrl", na.rm = TRUE)))

lib <- Matrix::colSums(X)
cat(sprintf("\n每细胞总 counts：中位数 %s，最小 %s\n",
            format(round(median(lib)), big.mark = ","),
            format(min(lib), big.mark = ",")))
cat(sprintf("基因表达率 >5%% 的基因数：%d / %d\n",
            sum(Matrix::rowMeans(X != 0) > 0.05), nrow(X)))

cat("\n输出文件：\n")
cat("  work/data/counts.mtx     (基因 x 细胞，原始 counts)\n")
cat("  work/data/genes.txt\n  work/data/barcodes.txt\n  work/data/metadata.csv\n")
