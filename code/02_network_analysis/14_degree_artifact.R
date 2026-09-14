# 检验"显著基因数"是否只是出度的函数（方法学诊断）
dir <- "outputs/results/vko_ms_microglia"

analyze <- function(tag) {
  s <- utils::read.csv(file.path(dir, paste0(tag, "_panel_summary.csv")), stringsAsFactors = FALSE)
  o <- utils::read.csv(file.path(dir, paste0(tag, "_outdegree.csv")), stringsAsFactors = FALSE)
  m <- merge(s[, c("ko_gene", "n_sig_FDR05", "is_negative_control")], o,
             by.x = "ko_gene", by.y = "gene")
  m <- m[order(m$outdegree), ]
  list(tag = tag, m = m)
}

for (tag in c("MSvC_Control", "MSvC_MS")) {
  a <- analyze(tag)
  m <- a$m
  cat(sprintf("\n===== %s =====\n", tag))
  cat(sprintf("%-9s %10s %12s %14s\n", "敲除基因", "出度", "显著基因数", "显著数-出度"))
  for (i in seq_len(nrow(m))) {
    cat(sprintf("%-9s %10d %12d %14d %s\n", m$ko_gene[i], m$outdegree[i], m$n_sig_FDR05[i],
                m$n_sig_FDR05[i] - m$outdegree[i],
                if (m$is_negative_control[i]) "[阴性对照]" else ""))
  }
  ## 只看出度 >0 且显著数 >1 的基因
  sub <- m[m$outdegree > 0 & m$n_sig_FDR05 > 1, ]
  if (nrow(sub) >= 3) {
    cat(sprintf("\n出度>0 且显著数>1 的基因：n=%d，n_sig 与 outdegree 的相关 rho = %.3f\n",
                nrow(sub), stats::cor(sub$outdegree, sub$n_sig_FDR05, method = "spearman")))
    cat(sprintf("差值 (n_sig - outdegree) 范围：%d 到 %d\n",
                min(sub$n_sig_FDR05 - sub$outdegree), max(sub$n_sig_FDR05 - sub$outdegree)))
  }
  zero <- m[m$outdegree == 0, ]
  cat(sprintf("\n出度为 0 的基因：%d 个，其显著基因数 = %s（空操作，统计量退化）\n",
              nrow(zero), paste(zero$n_sig_FDR05, collapse = ", ")))
}
