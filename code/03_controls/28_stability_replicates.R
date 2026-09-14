#!/usr/bin/env Rscript
# =============================================================================
# One stability replicate: rebuild the network from an independent cell
# subsample and run the full-gene virtual knockout.
#
# Usage:  Rscript 28_stability_replicates.R <tag> <rep> [n_cells] [n_genes]
#   tag      network tag, e.g. MSvC_Control or MSvC_MS
#   rep      replicate index (1, 2, 3, ...); controls the cell-subsampling seed
#   n_cells  cells to subsample per replicate (default 800)
#   n_genes  genes in the network (default 800)
#
# Output: work/stability/<tag>_rep<rep>.rds
#         list(network = WT, distances = D, outdegree = od, cells = <barcodes>)
#
# Each replicate is independent and single-threaded in practice, so replicates
# can be run in parallel (see run_stability_slurm.sh / run_stability_local.ps1).
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldNet) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: 28_stability_replicates.R <tag> <rep> [n_cells] [n_genes]")
tag     <- args[1]
rep_idx <- as.integer(args[2])
n_cells <- if (length(args) >= 3) as.integer(args[3]) else 800L
n_genes <- if (length(args) >= 4) as.integer(args[4]) else 800L

out_dir <- "work/stability"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("%s_rep%02d.rds", tag, rep_idx))
if (file.exists(out_file)) { cat("已存在，跳过：", out_file, "\n"); quit(status = 0) }

condition <- if (grepl("MS$", tag)) "MS" else "Control"
cat(sprintf("[%s rep %d] condition=%s, n_cells=%d, n_genes=%d\n",
            tag, rep_idx, condition, n_cells, n_genes))

## ---------------------------------------------------------------- input data
X <- Matrix::readMM("work/data/counts.mtx")
rownames(X) <- readLines("work/data/genes.txt")
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

## Same gene selection as the primary analysis: shared gene set (expression in
## >5% of all target cells) then the n_genes most variable, plus all panel genes.
PANEL <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA",
           "ACTA2","TSHR","PECAM1","YKT6")
cells_all <- intersect(meta$barcode[meta$celltype == "Microglia"], colnames(X))
Xall <- X[, cells_all, drop = FALSE]
rate <- Matrix::rowMeans(Xall != 0)
shared <- rownames(X)[rate > 0.05]
v <- as.numeric(Matrix::rowMeans(Xall * Xall) - Matrix::rowMeans(Xall)^2)
names(v) <- rownames(X)
rm(Xall); gc()
keep_genes <- union(names(sort(v[shared], decreasing = TRUE))[seq_len(n_genes)],
                    intersect(PANEL, shared))
rm(v); gc()

## Independent cell subsample for this replicate
set.seed(20260913 + rep_idx)
cells_cond <- intersect(meta$barcode[meta$condition == condition], colnames(X))
use <- if (length(cells_cond) > n_cells) sample(cells_cond, n_cells) else cells_cond

Xs <- X[intersect(keep_genes, rownames(X)), use, drop = FALSE]
Xs <- scTenifoldKnk::scQC(Xs, minLibSize = 500, removeOutlierCells = TRUE,
                          minPCT = 0, maxMTratio = 0.1, label = sprintf("%s_r%02d", tag, rep_idx))
Xs <- Xs[intersect(rownames(Xs), keep_genes), ]
drop <- grepl("^Rp[[:digit:]]|^Rpl|^Rps|^Mt-", rownames(Xs), ignore.case = TRUE)
Xs <- Xs[!drop, ]
Xs <- scTenifoldNet::cpmNormalization(Xs)
cat(sprintf("[%s rep %d] network input: %d genes x %d cells\n",
            tag, rep_idx, nrow(Xs), ncol(Xs)))

## ------------------------------------------------------------ build and knock out
t0 <- Sys.time()
set.seed(1)
W <- scTenifoldNet::makeNetworks(X = Xs, q = 0.90, nNet = 10, nCells = 500,
                                 scaleScores = TRUE, symmetric = FALSE, nComp = 3,
                                 nCores = 1)
set.seed(1)
W <- scTenifoldNet::tensorDecomposition(xList = W, K = 3, maxError = 1e-5,
                                        maxIter = 1000, nDecimal = 3)
W <- as.matrix(W$X)
diag(W) <- 0
W <- t(W)
cat(sprintf("[%s rep %d] network built in %.1f min (density %.1f%%)\n",
            tag, rep_idx, as.numeric(difftime(Sys.time(), t0, units = "mins")),
            100 * mean(W != 0)))

genes <- rownames(W)
D <- matrix(NA_real_, nrow = length(genes), ncol = ncol(W),
            dimnames = list(genes, colnames(W)))
t1 <- Sys.time()
for (i in seq_along(genes)) {
  g <- genes[i]
  KO <- W; KO[g, ] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(W, KO, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  D[g, DR$gene] <- DR$distance
}
cat(sprintf("[%s rep %d] knockouts done in %.1f min\n",
            tag, rep_idx, as.numeric(difftime(Sys.time(), t1, units = "mins"))))

saveRDS(list(network = W, distances = D,
             outdegree = Matrix::rowSums(W != 0),
             cells = colnames(Xs), genes = genes,
             params = list(tag = tag, rep = rep_idx, condition = condition,
                           n_cells = ncol(Xs), n_genes = nrow(Xs),
                           nNet = 10, nc_nCells = 500, q = 0.90, K = 3)),
        out_file)
cat("[", tag, " rep ", rep_idx, "] saved: ", out_file, "\n", sep = "")
