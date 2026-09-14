# 比较 Control 网络与 MS 网络中同一面板基因的扰动谱
# 目的：哪些基因的"敲除效应"在疾病中被重塑

suppressPackageStartupMessages({ library(Matrix) })

dir <- "outputs/results/vko_ms_microglia"
a <- as.matrix(utils::read.csv(file.path(dir, "MSvC_Control_distanceMatrix.csv"), row.names = 1, check.names = FALSE))
b <- as.matrix(utils::read.csv(file.path(dir, "MSvC_MS_distanceMatrix.csv"), row.names = 1, check.names = FALSE))

cat("Control 网络扰动矩阵：", dim(a)[1], "基因 x", dim(a)[2], "靶基因\n")
cat("MS 网络扰动矩阵：     ", dim(b)[1], "基因 x", dim(b)[2], "靶基因\n")

ko <- intersect(rownames(a), rownames(b))
genes <- intersect(colnames(a), colnames(b))
cat("共同敲除基因：", length(ko), " 共同被扰动基因：", length(genes), "\n\n")
a <- a[ko, genes, drop = FALSE]
b <- b[ko, genes, drop = FALSE]

neg <- c("ACTA2", "TSHR", "PECAM1", "YKT6")

res <- do.call(rbind, lapply(ko, function(g) {
  va <- a[g, ]; vb <- b[g, ]
  ok <- is.finite(va) & is.finite(vb)
  va <- va[ok]; vb <- vb[ok]
  rho <- suppressWarnings(stats::cor(va, vb, method = "spearman"))
  ta <- names(sort(va, decreasing = TRUE))[1:50]
  tb <- names(sort(vb, decreasing = TRUE))[1:50]
  jac <- length(intersect(ta, tb)) / length(union(ta, tb))
  ## 疾病网络中该基因的扰动强度（平均距离）
  data.frame(ko_gene = g,
             mean_dist_control = mean(va),
             mean_dist_MS = mean(vb),
             ratio_MS_vs_Control = mean(vb) / mean(va),
             rho_between_networks = rho,
             top50_jaccard = jac,
             is_negative_control = g %in% neg)
}))

res$rank_shift_score <- with(res, abs(scale(rank(mean_dist_MS)) - scale(rank(mean_dist_control))))
res <- res[order(-res$rho_between_networks), ]

cat("===== 每个敲除基因：两网络扰动谱的一致性 =====\n")
print(res, row.names = FALSE, digits = 3)

cat("\n===== 概览 =====\n")
cat(sprintf("rho 中位数：全部 %.2f ｜ 候选基因 %.2f ｜ 阴性对照 %.2f\n",
            median(res$rho_between_networks),
            median(res$rho_between_networks[!res$is_negative_control]),
            median(res$rho_between_networks[res$is_negative_control])))
cat(sprintf("疾病网络扰动强度／对照网络（中位比值）：%.2f\n", median(res$ratio_MS_vs_Control)))

utils::write.csv(res, file.path(dir, "network_comparison.csv"), row.names = FALSE)

png(file.path(dir, "network_comparison.png"), width = 1600, height = 900, res = 200)
op <- graphics::par(mar = c(10, 4, 3, 1))
bp <- graphics::barplot(res$rho_between_networks, names.arg = res$ko_gene, las = 2,
                        col = ifelse(res$is_negative_control, "grey70", "steelblue"),
                        main = "同一基因敲除效应在两套网络中的一致性（Spearman rho）",
                        ylab = "rho (Control vs MS 扰动谱)")
graphics::abline(h = median(res$rho_between_networks[res$is_negative_control]), lty = 2, col = "red")
graphics::par(op)
dev.off()

cat("\n已输出：network_comparison.csv / network_comparison.png\n")
