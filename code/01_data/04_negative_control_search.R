# 数据驱动：找表达率 5%-20%、跨细胞方差最低的基因，作为阴性对照候选
suppressPackageStartupMessages({ library(Matrix) })

X <- Matrix::readMM("work/data/counts.mtx")
genes <- readLines("work/data/genes.txt")
rownames(X) <- genes
colnames(X) <- readLines("work/data/barcodes.txt")

rate <- Matrix::rowMeans(X != 0)
keep <- which(rate >= 0.05 & rate <= 0.20)
cat("表达率 5%-20% 的基因数：", length(keep), "\n")

Y <- as.matrix(X[keep, ])
v <- apply(Y, 1, stats::var)
m <- rowMeans(Y)
ord <- order(v)

tab <- data.frame(
  gene = genes[keep][ord][1:25],
  expr_rate = sprintf("%.1f%%", 100 * rate[keep][ord][1:25]),
  mean_counts = round(m[ord][1:25], 2),
  variance = signif(v[ord][1:25], 2)
)
cat("\n===== 表达率 5%-20% 且方差最低的 25 个基因 =====\n")
print(tab, row.names = FALSE)
