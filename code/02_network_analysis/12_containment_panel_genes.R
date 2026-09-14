# 验证核心命题：显著基因集 ⊆ (被敲除基因的直接靶标 ∪ {自身})，且 |DR| ≤ 出度 + 1
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"

check <- function(tag) {
  res <- readRDS(file.path(dir, paste0(tag, "_raw_result.rds")))
  WT <- as.matrix(res$tensorNetworks$WT)
  DR <- res$diffRegulation                     # 主敲除（单基因模式）
  tab <- utils::read.csv(file.path(dir, paste0(tag, "_panel_allStats.csv")),
                         stringsAsFactors = FALSE)
  od <- Matrix::rowSums(WT != 0)

  cat(sprintf("\n===== %s =====\n", tag))
  cat(sprintf("%-9s %8s %8s %10s %12s %12s %s\n",
              "敲除基因", "出度", "显著数", "出度+1", "显著∈靶标∪自身", "靶标∩显著", "判定"))
  out <- do.call(rbind, lapply(unique(tab$ko_gene), function(g) {
    if (!g %in% rownames(WT)) return(NULL)
    sub <- tab[tab$ko_gene == g, ]
    sig <- sub$gene[sub$p.adj < 0.05]
    targets <- colnames(WT)[WT[g, ] != 0]
    allowed <- union(targets, g)
    inside <- all(sig %in% allowed)
    n_ov <- length(intersect(sig, targets))
    data.frame(ko = g, od = od[g], n_sig = length(sig), bound = od[g] + 1,
               all_inside = inside, n_target_sig = n_ov,
               frac_targets_sig = if (length(targets) > 0) round(100 * n_ov / length(targets), 1) else NA)
  }))
  out <- out[order(out$od), ]
  for (i in seq_len(nrow(out))) {
    cat(sprintf("%-9s %8d %8d %10d %12s %12d (%.1f%%) %s\n",
                out$ko[i], out$od[i], out$n_sig[i], out$bound[i],
                ifelse(out$all_inside[i], "是", "**否**"),
                out$n_target_sig[i], out$frac_targets_sig[i],
                ifelse(out$n_sig[i] <= out$bound[i], "OK", "**越界**")))
  }
  cat(sprintf("\n汇总：%d/%d 个基因的显著集完全落在'靶标∪自身'内；%d/%d 满足 |DR| ≤ 出度+1\n",
              sum(out$all_inside), nrow(out), sum(out$n_sig <= out$bound), nrow(out)))
  out
}

a <- check("MSvC_Control")
b <- check("MSvC_MS")
write.csv(a, file.path(dir, "dr_bound_Control.csv"), row.names = FALSE)
write.csv(b, file.path(dir, "dr_bound_MS.csv"), row.names = FALSE)
