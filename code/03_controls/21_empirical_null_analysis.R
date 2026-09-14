# 基于全基因经验零分布，对面板基因做经验 p 值检验
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"
panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")
neg <- c("ACTA2","TSHR","PECAM1","YKT6")

analyze <- function(tag) {
  D <- readRDS(file.path(dir, paste0(tag, "_fullKO_distances.rds")))
  od <- read.csv(file.path(dir, paste0(tag, "_outdegree.csv")), stringsAsFactors = FALSE)
  odv <- setNames(od$outdegree, od$gene)[rownames(D)]

  ## 特异性 z：对每个靶基因，跨全部 800 次敲除标准化
  mu <- colMeans(D, na.rm = TRUE)
  s <- apply(D, 2, stats::sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- NA
  Z <- sweep(sweep(D, 2, mu, "-"), 2, s, "/")
  Z[!is.finite(Z)] <- NA
  diag(Z) <- NA

  maxz <- apply(Z, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  maxz[!is.finite(maxz)] <- NA

  ## 经验零分布：仅用非面板基因
  null_genes <- setdiff(rownames(D), c(panel, neg))
  null_maxz <- maxz[null_genes]
  null_od <- odv[null_genes]

  eval_gene <- function(g) {
    mz <- unname(maxz[g]); o <- unname(odv[g])
    if (is.na(mz)) return(c(NA, NA, NA, NA))
    p_all <- (1 + sum(null_maxz >= mz, na.rm = TRUE)) / (1 + sum(!is.na(null_maxz)))
    ## 出度匹配（±20%）
    lo <- o * 0.8; hi <- o * 1.2
    sel <- !is.na(null_maxz) & !is.na(null_od) & null_od >= lo & null_od <= hi
    p_matched <- if (sum(sel) >= 20) (1 + sum(null_maxz[sel] >= mz)) / (1 + sum(sel)) else NA
    c(maxz = round(mz, 2), p_all = signif(p_all, 3),
      p_matched = signif(p_matched, 3), n_matched = sum(sel))
  }

  rows <- lapply(c(panel, neg), eval_gene)
  res <- data.frame(
    ko_gene = c(panel, neg),
    role = c(rep("candidate", length(panel)), rep("negative_control", length(neg))),
    outdegree = as.integer(odv[c(panel, neg)]),
    maxz = as.numeric(vapply(rows, function(x) x[["maxz"]], numeric(1))),
    p_all = as.numeric(vapply(rows, function(x) x[["p_all"]], numeric(1))),
    p_matched = as.numeric(vapply(rows, function(x) x[["p_matched"]], numeric(1))),
    n_matched = as.integer(vapply(rows, function(x) x[["n_matched"]], numeric(1))),
    stringsAsFactors = FALSE
  )
  res$p_all_BH <- stats::p.adjust(res$p_all, method = "fdr")
  res <- res[order(-res$maxz), ]
  list(res = res, null_maxz = null_maxz, maxz = maxz)
}

out <- list()
for (tag in c("MSvC_Control", "MSvC_MS")) {
  a <- analyze(tag)
  out[[tag]] <- a
  cat(sprintf("\n===== %s：经验零分布检验 =====\n", tag))
  cat(sprintf("零分布（%d 个非面板基因）max_z：中位 %.2f，95 分位 %.2f，最大 %.2f\n",
              length(a$null_maxz), median(a$null_maxz, na.rm = TRUE),
              quantile(a$null_maxz, 0.95, na.rm = TRUE), max(a$null_maxz, na.rm = TRUE)))
  print(a$res, row.names = FALSE)
  write.csv(a$res, file.path(dir, paste0("empirical_pvalue_", sub("MSvC_", "", tag), ".csv")),
            row.names = FALSE)
}

## 画图
png(file.path(dir, "empirical_null_plot.png"), width = 1700, height = 900, res = 180)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (tag in c("MSvC_Control", "MSvC_MS")) {
  a <- out[[tag]]
  h <- graphics::hist(a$null_maxz, breaks = 40, plot = FALSE)
  graphics::plot(h, col = "grey85", border = "white",
                 main = paste0(sub("MSvC_", "", tag), "：经验零分布"),
                 xlab = "max z（非面板基因）")
  sig <- a$res[a$res$role == "candidate" & a$res$maxz > 2, ]
  graphics::abline(v = quantile(a$null_maxz, 0.95, na.rm = TRUE), lty = 2, col = "red")
  if (nrow(sig) > 0) {
    graphics::points(sig$maxz, rep(0, nrow(sig)), pch = 25, bg = "steelblue", cex = 1.4)
    graphics::text(sig$maxz, 0, labels = sig$ko_gene, pos = 3, cex = 0.7, srt = 45)
  }
}
graphics::par(op); dev.off()
cat("\n完成：empirical_pvalue_*.csv / empirical_null_plot.png\n")
