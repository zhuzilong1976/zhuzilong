# 诊断：面板基因在推断网络中的出度（出度为 0 = 虚拟敲除是空操作）
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })

rds <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(rds)) rds <- "work/test_output_real/MSmg_Control_smoke_raw_result.rds"
res <- readRDS(rds)
WT <- res$tensorNetworks$WT

cat(sprintf("网络规模：%d x %d，非零边 %s（密度 %.2f%%）\n",
            nrow(WT), ncol(WT), format(length(WT@x), big.mark = ","),
            100 * length(WT@x) / (nrow(WT) * ncol(WT))))

outdeg <- Matrix::rowSums(WT != 0)
outdeg <- sort(outdeg, decreasing = TRUE)
cat(sprintf("出度分布：中位 %.0f，最大 %d，出度为 0 的基因 %d 个（%.1f%%）\n",
            median(outdeg), max(outdeg), sum(outdeg == 0), 100 * mean(outdeg == 0)))

panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA",
           "ACTA2","TSHR","PECAM1","YKT6")

cat("\n===== 面板基因的出度与排名 =====\n")
cat(sprintf("%-9s %9s %9s %s\n", "基因", "出度", "百分位", "评估"))
for (g in panel) {
  if (!g %in% rownames(WT)) { cat(sprintf("%-9s 不在网络中\n", g)); next }
  d <- outdeg[g]
  pct <- round(100 * mean(outdeg <= d))
  verdict <- if (d == 0) "★ 出度为 0：敲除无效！" else if (pct < 50) "低（扰动可能很弱）" else "正常"
  cat(sprintf("%-9s %9d %8d%% %s\n", g, d, pct, verdict))
}
