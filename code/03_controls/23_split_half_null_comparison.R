# 技术对照：用对照组两个半样本网络，测量"无疾病差异时"能产生多少假阳性
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })
dir <- "outputs/results/vko_ms_microglia"
ndir <- "outputs/results/technical_null"

run_comp <- function(A, B, label) {
  g <- intersect(rownames(A), rownames(B))
  A <- A[g, g]; B <- B[g, g]
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(A, B, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  d2 <- DR$distance^2
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("基因数 %d；p<0.05 的 %d 个（%.1f%%）；FDR<0.05 的 %d 个（%.1f%%）\n",
              nrow(DR), sum(DR$p.value < 0.05), 100 * mean(DR$p.value < 0.05),
              sum(DR$p.adj < 0.05), 100 * mean(DR$p.adj < 0.05)))
  cat(sprintf("距离中位 %.2e，最大 %.2e；最大基因占信号量 %.1f%%\n",
              median(DR$distance), max(DR$distance), 100 * max(d2) / sum(d2)))
  cat("Top 10：", paste(head(DR$gene[order(DR$p.value)], 10), collapse = ", "), "\n")
  list(DR = DR, n_sig = sum(DR$p.value < 0.05), frac_sig = mean(DR$p.value < 0.05),
       top = head(DR$gene[order(DR$p.value)], 100))
}

Aa <- readRDS(file.path(ndir, "CtrlHalfA_network.rds"))
Bb <- readRDS(file.path(ndir, "CtrlHalfB_network.rds"))
null_res <- run_comp(Aa, Bb, "技术对照：对照组半样本 A vs B（无疾病差异）")

obs <- readRDS(file.path(dir, "condition_comparison.rds"))
obs_res <- run_comp(obs$A, obs$B, "实际比较：Control vs MS 网络")

cat("\n\n================= 判定 =================\n")
cat(sprintf("技术对照假阳性率：%.1f%%（预期约 5%%）\n", 100 * null_res$frac_sig))
cat(sprintf("实际比较显著率：  %.1f%%\n", 100 * obs_res$frac_sig))
cat(sprintf("倍数：%.1f×\n", obs_res$frac_sig / null_res$frac_sig))
ov <- length(intersect(null_res$top, obs_res$top)) / length(union(null_res$top, obs_res$top))
cat(sprintf("Top100 重叠 Jaccard：%.3f（越接近 0 说明两组结果越独立）\n", ov))
cat("\n判读标准：\n")
cat("  若实际显著率 ≈ 对照假阳性率 → 疾病 vs 对照的网络差异无法与测量噪声区分\n")
cat("  若实际显著率明显更高且 Top 基因不同 → 存在真实的疾病相关网络重排\n")

png(file.path(dir, "null_vs_observed.png"), width = 1500, height = 800, res = 180)
op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
graphics::hist(null_res$DR$distance, breaks = 60, col = "grey75", border = "white",
               main = "技术对照（同条件两个半样本）", xlab = "distance")
graphics::hist(obs_res$DR$distance, breaks = 60, col = "firebrick", border = "white",
               main = "实际比较（Control vs MS）", xlab = "distance")
graphics::par(op); dev.off()
cat("\n完成：null_vs_observed.png\n")
