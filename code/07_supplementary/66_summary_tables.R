# B 方案最终产出：特异性指标 + 高信号基因的靶标列表 + 跨条件比较图
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"

read_mat <- function(tag) as.matrix(utils::read.csv(file.path(dir, paste0(tag, "_distanceMatrix.csv")),
                                                    row.names = 1, check.names = FALSE))
ctrl <- read_mat("MSvC_Control"); ms <- read_mat("MSvC_MS")
mc <- utils::read.csv(file.path(dir, "specificity_metric_Control.csv"), stringsAsFactors = FALSE)
mm <- utils::read.csv(file.path(dir, "specificity_metric_MS.csv"), stringsAsFactors = FALSE)

spec_z <- function(mat) {
  M <- mat
  mu <- colMeans(M); s <- apply(M, 2, stats::sd); s[s == 0] <- NA
  Z <- sweep(sweep(M, 2, mu, "-"), 2, s, "/"); Z[!is.finite(Z)] <- 0; Z
}
Zc <- spec_z(ctrl); Zm <- spec_z(ms)

## ---- 1. 跨条件比较 ----
cmp <- merge(mc[, c("ko_gene", "outdegree", "max_z", "n_z2", "is_negative_control")],
             mm[, c("ko_gene", "outdegree", "max_z", "n_z2")],
             by = "ko_gene", suffixes = c("_Control", "_MS"))
cmp <- cmp[order(-pmax(cmp$max_z_Control, cmp$max_z_MS, na.rm = TRUE)), ]
cmp$delta_max_z <- round(cmp$max_z_MS - cmp$max_z_Control, 2)
cat("===== 跨条件比较（按较强的一侧排序）=====\n")
print(cmp, row.names = FALSE)
utils::write.csv(cmp, file.path(dir, "specificity_metric_comparison.csv"), row.names = FALSE)

## ---- 2. 高信号基因的靶标（z>2）----
top_targets <- function(Z, g, k = 15) {
  v <- Z[g, ]; v <- v[names(v) != g]
  v <- sort(v, decreasing = TRUE)
  v[seq_len(min(k, length(v)))]
}
cat("\n\n===== 高信号基因的特异性靶标（Top 15，z 值）=====\n")
for (nm in c("Control", "MS")) {
  Z <- if (nm == "Control") Zc else Zm
  m <- if (nm == "Control") mc else mm
  strong <- m$ko_gene[!is.na(m$max_z) & m$max_z > 2 & !m$is_negative_control]
  cat(sprintf("\n----- %s 网络：%d 个基因有信号 -----\n", nm, length(strong)))
  for (g in strong) {
    tt <- top_targets(Z, g)
    cat(sprintf("\n[%s] max_z=%.2f, 出度=%d\n   %s\n", g, max(tt), m$outdegree[m$ko_gene == g],
                paste(sprintf("%s(%.1f)%s", names(tt), tt, ifelse(tt > 2, "*", "")), collapse = "  ")))
  }
}

## ---- 3. 比较图 ----
png(file.path(dir, "specificity_metric_plot.png"), width = 1700, height = 900, res = 180)
op <- graphics::par(mar = c(9, 4, 3, 1))
M <- rbind(Control = mc$max_z[match(cmp$ko_gene, mc$ko_gene)],
           MS = mm$max_z[match(cmp$ko_gene, mm$ko_gene)])
M[is.na(M)] <- 0
graphics::barplot(M, beside = TRUE, names.arg = cmp$ko_gene, las = 2,
                  col = c("steelblue", "firebrick"),
                  main = "特异性扰动强度（max z）：MS vs Control 网络",
                  ylab = "max z", legend.text = c("Control", "MS"),
                  args.legend = list(x = "topright", bty = "n"))
graphics::abline(h = 2, lty = 2, col = "grey40")
graphics::par(op)
dev.off()

cat("\n\n完成：specificity_metric_comparison.csv / specificity_metric_plot.png\n")
