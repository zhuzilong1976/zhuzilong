# 构建单个条件网络（复用与主分析完全一致的流程），用于技术对照（split-half null）
# 用法：Rscript work/build_network.R <输出标签> <条件> <A|B|all>
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldNet) })

args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]; cond <- args[2]; half <- args[3]
dir <- "outputs/results/technical_null"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

X <- Matrix::readMM("work/data/counts.mtx")
rownames(X) <- readLines("work/data/genes.txt")
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

## ---- 与主分析一致的基因集与质控 ----
cells_target <- intersect(meta$barcode[meta$celltype == "Microglia"], colnames(X))
rate_all <- Matrix::rowMeans(X[, cells_target] != 0)
shared <- rownames(X)[rate_all > 0.05]
Xt <- X[, cells_target, drop = FALSE]
v <- as.numeric(Matrix::rowMeans(Xt * Xt) - Matrix::rowMeans(Xt)^2)
names(v) <- rownames(X)
rm(Xt); gc()
top <- names(sort(v[shared], decreasing = TRUE))[1:800]
panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA",
           "ACTA2","TSHR","PECAM1","YKT6")
kg <- union(top, intersect(panel, shared))
rm(v); gc()

cells <- intersect(meta$barcode[meta$condition == cond], colnames(X))
set.seed(20260913)
ord <- sample(cells)
halfA <- ord[seq(1, length(ord), by = 2)]
halfB <- ord[seq(2, length(ord), by = 2)]
use <- switch(half, A = halfA, B = halfB, all = cells)
use <- intersect(use, colnames(X))

Xs <- X[intersect(kg, rownames(X)), use, drop = FALSE]
Xs <- scTenifoldKnk::scQC(Xs, minLibSize = 500, removeOutlierCells = TRUE,
                          minPCT = 0, maxMTratio = 0.1, label = tag)
Xs <- Xs[intersect(rownames(Xs), kg), ]
keep <- !grepl("^Rp[[:digit:]]|^Rpl|^Rps|^Mt-", rownames(Xs), ignore.case = TRUE)
Xs <- Xs[keep, ]
Xs <- scTenifoldNet::cpmNormalization(Xs)
cat(sprintf("[%s] 输入 %d 基因 x %d 细胞\n", tag, nrow(Xs), ncol(Xs)))

set.seed(1)
W <- scTenifoldNet::makeNetworks(X = Xs, q = 0.9, nNet = 10, nCells = 500,
                                 scaleScores = TRUE, symmetric = FALSE, nComp = 3,
                                 nCores = max(1, parallel::detectCores() - 1))
set.seed(1)
W <- scTenifoldNet::tensorDecomposition(xList = W, K = 3, maxError = 1e-5, maxIter = 1000,
                                        nDecimal = 3)
W <- as.matrix(W$X)
diag(W) <- 0
W <- t(W)
saveRDS(W, file.path(dir, paste0(tag, "_network.rds")))
cat(sprintf("[%s] 网络完成并保存：%d x %d\n", tag, nrow(W), ncol(W)))
