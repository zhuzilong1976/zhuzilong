# 全基因经验零分布：复用已建网络，对全部网络基因做虚拟敲除
# 用法：Rscript work/empirical_null.R MSvC_Control
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })

tag <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(tag)) tag <- "MSvC_Control"
dir <- "outputs/results/vko_ms_microglia"

res <- readRDS(file.path(dir, paste0(tag, "_raw_result.rds")))
WT <- as.matrix(res$tensorNetworks$WT)
genes <- rownames(WT)
cat(sprintf("[%s] 网络 %d x %d，非零 %s\n", tag, nrow(WT), ncol(WT),
            format(sum(WT != 0), big.mark = ",")))

od <- Matrix::rowSums(WT != 0)
D <- matrix(NA_real_, nrow = length(genes), ncol = ncol(WT),
            dimnames = list(genes, colnames(WT)))

t0 <- Sys.time()
for (i in seq_along(genes)) {
  g <- genes[i]
  KO <- WT
  KO[g, ] <- 0
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  D[g, DR$gene] <- DR$distance
  if (i %% 100 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %d/%d 完成，已用 %.1f 分钟，预计剩余 %.1f 分钟\n",
                i, length(genes), el, el / i * (length(genes) - i)))
  }
}
cat(sprintf("[%s] 全部完成，用时 %.1f 分钟\n", tag,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

saveRDS(D, file.path(dir, paste0(tag, "_fullKO_distances.rds")))
cat("已保存：", paste0(tag, "_fullKO_distances.rds"), "\n")
