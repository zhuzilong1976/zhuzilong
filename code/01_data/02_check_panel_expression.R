# 在真实数据里核对基因面板，并挑选阴性对照基因

suppressPackageStartupMessages({ library(Matrix) })

X <- Matrix::readMM("work/data/counts.mtx")
genes <- readLines("work/data/genes.txt")
barcodes <- readLines("work/data/barcodes.txt")
rownames(X) <- genes
colnames(X) <- barcodes
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")

expr_rate <- Matrix::rowMeans(X != 0)
mean_expr <- Matrix::rowMeans(X)

cat("===== 候选基因面板在 MS 小胶质中的表达 =====\n")
cat(sprintf("%-10s %10s %12s %10s\n", "基因", "表达率", "平均counts", "是否在矩阵"))
for (g in panel) {
  if (g %in% genes) {
    cat(sprintf("%-10s %9.1f%% %12.2f %10s\n", g, 100 * expr_rate[g], mean_expr[g], "是"))
  } else {
    cat(sprintf("%-10s %10s %12s %10s\n", g, "-", "-", "否 <== 缺失"))
  }
}

## 对照条件与疾病条件的表达率对比（看基因是否随疾病变化）
ms_cells <- meta$barcode[meta$condition == "MS"]
ct_cells <- meta$barcode[meta$condition == "Control"]
Xms <- X[, ms_cells]; Xct <- X[, ct_cells]
rate_ms <- Matrix::rowMeans(Xms != 0)
rate_ct <- Matrix::rowMeans(Xct != 0)

cat("\n===== 候选基因：MS vs Control 表达率 =====\n")
cat(sprintf("%-10s %10s %10s %10s\n", "基因", "MS", "Control", "差值"))
for (g in panel) {
  if (g %in% genes) {
    cat(sprintf("%-10s %9.1f%% %9.1f%% %+9.1f%%\n", g,
                100 * rate_ms[g], 100 * rate_ct[g], 100 * (rate_ms[g] - rate_ct[g])))
  }
}

## ---- 阴性对照候选：在 1%–5% 细胞中表达、且与神经免疫无关的谱系基因 ----
candidates <- c("ALB","HBB","HBA1","HBA2","INS","GCG","TSHR","KRT5","KRT14","MYH7",
                "CD3E","CD19","MS4A1","NKG7","EPCAM","PTPRC")
cat("\n===== 阴性对照候选（应低表达/不表达）=====\n")
cat(sprintf("%-10s %10s %10s %10s\n", "基因", "MS", "Control", "可用性"))
for (g in candidates) {
  if (g %in% genes) {
    note <- if (rate_ms[g] < 0.05 && rate_ct[g] < 0.05) "<== 可用作阴性对照" else if (rate_ms[g] < 0.1) "偏高" else "不适合"
    cat(sprintf("%-10s %9.1f%% %9.1f%%   %s\n", g, 100 * rate_ms[g], 100 * rate_ct[g], note))
  } else {
    cat(sprintf("%-10s %10s %10s   不在矩阵中\n", g, "-", "-"))
  }
}

## ---- 客观筛选：表达率 1%–5% 且方差最低的基因 ----
ok <- which(expr_rate > 0.01 & expr_rate < 0.05)
v <- apply(as.matrix(X[ok, ]), 1, stats::var)
ord <- order(v)
cat("\n===== 数据驱动的阴性对照候选（表达率 1%-5%，方差最低的 15 个）=====\n")
print(data.frame(
  gene = genes[ok][ord[1:15]],
  expr_rate = sprintf("%.1f%%", 100 * expr_rate[ok][ord[1:15]]),
  mean_counts = round(mean_expr[ok][ord[1:15]], 3)
), row.names = FALSE)
