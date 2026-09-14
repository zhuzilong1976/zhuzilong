# D 方案：直接比较"对照网络 vs 疾病网络"（不依赖单基因虚拟敲除）
# 复用已建好的两个条件网络，参数完全一致

suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })
dir <- "outputs/results/vko_ms_microglia"

rc <- readRDS(file.path(dir, "MSvC_Control_raw_result.rds"))
rm_ <- readRDS(file.path(dir, "MSvC_MS_raw_result.rds"))
A <- as.matrix(rc$tensorNetworks$WT)   # 对照网络
B <- as.matrix(rm_$tensorNetworks$WT)  # 疾病网络
cat("对照网络：", nrow(A), "x", ncol(A), " 非零", format(sum(A != 0), big.mark = ","), "\n")
cat("疾病网络：", nrow(B), "x", ncol(B), " 非零", format(sum(B != 0), big.mark = ","), "\n")

g <- intersect(rownames(A), rownames(B))
cat("共同基因：", length(g), "\n")
A <- A[g, g]; B <- B[g, g]

## 两个网络本身的相关性（整体相似度）
cat(sprintf("两网络整体 Spearman 相关：%.3f\n", stats::cor(as.vector(A), as.vector(B), method = "spearman")))
cat(sprintf("平均绝对差异：%.2e（对照均值 %.2e）\n",
            mean(abs(A - B)), mean(abs(A))))

## ---- 两条件比较 ----
set.seed(1)
MA <- scTenifoldNet::manifoldAlignment(A, B, d = 2)
DR <- scTenifoldKnk::dRegulation(MA)
DR <- DR[order(DR$p.value), ]

cat("\n===== 差异调控总体情况 =====\n")
cat(sprintf("基因数 %d；p<0.05 的 %d 个（%.1f%%）；FDR<0.05 的 %d 个（%.1f%%）\n",
            nrow(DR), sum(DR$p.value < 0.05), 100 * mean(DR$p.value < 0.05),
            sum(DR$p.adj < 0.05), 100 * mean(DR$p.adj < 0.05)))
d2 <- DR$distance^2
cat(sprintf("距离：中位 %.2e，最大 %.2e；最大基因占信号量 %.1f%%\n",
            median(DR$distance), max(DR$distance), 100 * max(d2) / sum(d2)))

cat("\nTop 20 差异调控基因：\n")
print(head(DR[, c("gene", "distance", "Z", "FC", "p.value", "p.adj")], 20), row.names = FALSE, digits = 3)

## ---- 面板基因在全部基因中的位置 ----
panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")
neg <- c("ACTA2","TSHR","PECAM1","YKT6")
DR$rank <- seq_len(nrow(DR))
DR$pct <- 100 * DR$rank / nrow(DR)
cat("\n===== 面板基因的排名与 p 值 =====\n")
sel <- DR[DR$gene %in% c(panel, neg), c("gene", "distance", "p.value", "p.adj", "rank", "pct")]
sel$role <- ifelse(sel$gene %in% neg, "阴性对照", "候选基因")
print(sel[order(sel$rank), ], row.names = FALSE, digits = 3)

cat("\n===== 面板基因 vs 其他基因的整体比较 =====\n")
inpanel <- DR$gene %in% panel
cat(sprintf("候选基因 p 值中位数 %.3f ｜ 其余基因 %.3f（Wilcoxon p = %.3g）\n",
            median(DR$p.value[inpanel]), median(DR$p.value[!inpanel]),
            suppressWarnings(stats::wilcox.test(DR$p.value[inpanel], DR$p.value[!inpanel])$p.value)))
cat(sprintf("候选基因平均排名百分位 %.1f%% ｜ 其余 %.1f%%\n",
            mean(DR$pct[inpanel]), mean(DR$pct[!inpanel])))

utils::write.csv(DR, file.path(dir, "condition_comparison_DR.csv"), row.names = FALSE)
saveRDS(list(A = A, B = B, MA = MA, DR = DR), file.path(dir, "condition_comparison.rds"))

## ---- 图 ----
png(file.path(dir, "condition_comparison.png"), width = 1700, height = 800, res = 180)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
graphics::hist(log10(DR$distance + 1e-12), breaks = 60, col = "steelblue", border = "white",
               main = "差异调控距离分布（log10）", xlab = "log10(distance)")
graphics::abline(v = log10(min(DR$distance[DR$p.adj < 0.05]) + 1e-12), lty = 2, col = "red")
qq <- stats::qqnorm(DR$Z, plot.it = FALSE)
graphics::plot(qq$x, qq$y, pch = 16, cex = 0.5,
               col = ifelse(DR$p.adj < 0.05, "red", "grey50"),
               main = "Z 值 QQ 图", xlab = "理论分位数", ylab = "观测 Z")
graphics::abline(0, 1, lty = 2)
graphics::par(op); dev.off()

cat("\n完成：condition_comparison_DR.csv / condition_comparison.png / condition_comparison.rds\n")
