# 参数敏感性：检验"DR ⊆ {KO} ∪ 直接靶标"是否依赖 nc_lambda / nc_q
# 用法：Rscript work/param_sensitivity.R <tag> <lambda> <q> [基因数] [每网络细胞数]
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldNet) })

args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]
lambda <- as.numeric(args[2])
qq <- as.numeric(args[3])
nGenes <- if (length(args) >= 4) as.integer(args[4]) else 400L
nCells <- if (length(args) >= 5) as.integer(args[5]) else 400L
nNet <- 10L
outdir <- "outputs/results/param_sensitivity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[%s] lambda=%.2f q=%.2f genes=%d cells/net=%d nets=%d\n",
            tag, lambda, qq, nGenes, nCells, nNet))

## scTenifoldKnk 内部的 strictDirection（未导出，此处按源码复制）
strictDirection <- function(X, lambda = 1) {
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  O <- (((1 - lambda) * X) + (lambda * S))
  Matrix::Matrix(O)
}

## ---- 数据与质控（与主分析一致）----
X <- Matrix::readMM("work/data/counts.mtx")
rownames(X) <- readLines("work/data/genes.txt")
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

cells_target <- intersect(meta$barcode[meta$celltype == "Microglia"], colnames(X))
Xt <- X[, cells_target, drop = FALSE]
rate_all <- Matrix::rowMeans(Xt != 0)
shared <- rownames(X)[rate_all > 0.05]
v <- as.numeric(Matrix::rowMeans(Xt * Xt) - Matrix::rowMeans(Xt)^2)
names(v) <- rownames(X)
rm(Xt); gc()
top <- names(sort(v[shared], decreasing = TRUE))[1:nGenes]
panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")
kg <- union(top, intersect(panel, shared))
rm(v); gc()

use <- intersect(meta$barcode[meta$condition == "Control"], colnames(X))
Xs <- X[intersect(kg, rownames(X)), use, drop = FALSE]
Xs <- scTenifoldKnk::scQC(Xs, minLibSize = 500, removeOutlierCells = TRUE,
                          minPCT = 0, maxMTratio = 0.1, label = tag)
Xs <- Xs[intersect(rownames(Xs), kg), ]
keep <- !grepl("^Rp[[:digit:]]|^Rpl|^Rps|^Mt-", rownames(Xs), ignore.case = TRUE)
Xs <- Xs[keep, ]
Xs <- scTenifoldNet::cpmNormalization(Xs)
cat(sprintf("[%s] 输入 %d 基因 x %d 细胞\n", tag, nrow(Xs), ncol(Xs)))

## ---- 建网（含参数）----
t0 <- Sys.time()
set.seed(1)
W <- makeNetworks(X = Xs, q = qq, nNet = nNet, nCells = nCells, scaleScores = TRUE,
                  symmetric = FALSE, nComp = 3, nCores = max(1, parallel::detectCores() - 1))
set.seed(1)
W <- tensorDecomposition(xList = W, K = 3, maxError = 1e-5, maxIter = 1000, nDecimal = 3)
W <- W$X
W <- strictDirection(W, lambda = lambda)
W <- as.matrix(W)
diag(W) <- 0
W <- t(W)
cat(sprintf("[%s] 建网完成，用时 %.1f 分钟；密度 %.1f%%\n", tag,
            as.numeric(difftime(Sys.time(), t0, units = "mins")), 100 * mean(W != 0)))

## ---- 全基因虚拟敲除 ----
genes <- rownames(W)
D <- matrix(NA_real_, nrow = length(genes), ncol = ncol(W), dimnames = list(genes, colnames(W)))
t1 <- Sys.time()
for (i in seq_along(genes)) {
  g <- genes[i]
  KO <- W; KO[g, ] <- 0
  set.seed(1)
  MA <- manifoldAlignment(W, KO, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  D[g, DR$gene] <- DR$distance
}
cat(sprintf("[%s] 敲除完成，用时 %.1f 分钟\n", tag,
            as.numeric(difftime(Sys.time(), t1, units = "mins"))))

## ---- 包含关系检验 ----
od <- Matrix::rowSums(W != 0)
res <- do.call(rbind, lapply(genes, function(g) {
  d <- D[g, ]
  FC <- d^2 / mean(d^2, na.rm = TRUE)
  p <- stats::pchisq(FC, df = 1, lower.tail = FALSE)
  padj <- stats::p.adjust(p, method = "fdr")
  sig <- names(padj)[padj < 0.05]
  targets <- colnames(W)[W[g, ] != 0]
  data.frame(ko = g, outdegree = od[g], n_sig = length(sig),
             all_inside = all(sig %in% union(targets, g)),
             frac_sig = mean(padj < 0.05))
}))
res$lambda <- lambda; res$q <- qq
utils::write.csv(res, file.path(outdir, paste0(tag, "_containment.csv")), row.names = FALSE)
saveRDS(list(net = W, D = D, res = res), file.path(outdir, paste0(tag, ".rds")))

cat(sprintf("\n===== %s 结果 =====\n", tag))
cat(sprintf("参数：lambda=%.2f  q=%.2f\n", lambda, qq))
cat(sprintf("网络密度：%.1f%%；出度中位 %d（最小 %d，最大 %d）\n",
            100 * mean(W != 0), median(od), min(od), max(od)))
cat(sprintf("包含关系合规：%d / %d = %.1f%%\n", sum(res$all_inside), nrow(res),
            100 * mean(res$all_inside)))
cat(sprintf("|DR| ≤ 出度+1：%d / %d = %.1f%%\n", sum(res$n_sig <= res$outdegree + 1),
            nrow(res), 100 * mean(res$n_sig <= res$outdegree + 1)))
cat(sprintf("显著基因数中位 %d；整体显著率 %.2f%%\n", median(res$n_sig), 100 * mean(res$frac_sig)))
bad <- res[!res$all_inside, ]
if (nrow(bad) > 0) cat(sprintf("违例 %d 个，其中出度为 0 的 %d 个\n", nrow(bad), sum(bad$outdegree == 0)))
