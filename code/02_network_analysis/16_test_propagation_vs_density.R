# 决定性实验：网络稀疏度是否决定虚拟敲除的传播能力
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })

res <- readRDS("outputs/results/vko_ms_microglia/MSvC_Control_raw_result.rds")
WT0 <- as.matrix(res$tensorNetworks$WT)
cat(sprintf("网络：%d x %d；非零 %s（密度 %.1f%%）\n",
            nrow(WT0), ncol(WT0), format(sum(WT0 != 0), big.mark = ","),
            100 * mean(WT0 != 0)))
cat("|权重| 分位数：", paste(round(quantile(abs(WT0), c(.5, .9, .95, .99, .999)), 5), collapse = " | "), "\n\n")

target <- "CSF1R"
cat(sprintf("%-10s %10s %12s %10s %12s\n", "保留比例", "非零边", "显著基因数", "KO自身FC", "KO占信号比"))

for (q in c(0, 0.5, 0.9, 0.95, 0.99)) {
  W <- WT0
  if (q > 0) {
    thr <- stats::quantile(abs(W), q)
    W[abs(W) < thr] <- 0
  }
  KO <- W
  KO[target, ] <- 0
  MA <- try(scTenifoldNet::manifoldAlignment(W, KO, d = 2), silent = TRUE)
  if (inherits(MA, "try-error")) { cat(sprintf("%-10s 失败\n", q)); next }
  DR <- scTenifoldKnk::dRegulation(MA)
  d2 <- DR$distance^2
  ko_fc <- d2[DR$gene == target] / mean(d2)
  cat(sprintf("top %-6.0f%% %10s %12d %10.1f %11.1f%%\n",
              100 * (1 - q), format(sum(W != 0), big.mark = ","),
              sum(DR$p.adj < 0.05), ko_fc, 100 * max(d2) / sum(d2)))
}

cat("\n===== 对照：敲除一个与网络无关的基因 =====\n")
W <- WT0
thr <- stats::quantile(abs(W), 0.95); W[abs(W) < thr] <- 0
for (g in c("CSF1R", "TREM2")) {
  KO <- W; KO[g, ] <- 0
  MA <- scTenifoldNet::manifoldAlignment(W, KO, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  d2 <- DR$distance^2
  cat(sprintf("  %-8s 显著基因数 %3d，KO自身 FC %8.1f，KO占信号 %.1f%%\n",
              g, sum(DR$p.adj < 0.05), d2[DR$gene == g] / mean(d2), 100 * max(d2) / sum(d2)))
}
