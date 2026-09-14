# 诊断：虚拟敲除的效应是否在网络中传播
dir <- "outputs/results/vko_ms_microglia"
st <- utils::read.csv(file.path(dir, "MSvC_Control_panel_allStats.csv"), stringsAsFactors = FALSE)

cat("===== 每个敲除基因：传播能力诊断 =====\n")
cat(sprintf("%-9s %7s %12s %12s %10s %10s\n",
            "敲除基因", "出度档", "最大FC", "中位FC", "FC>3.84数", "FDR<0.05"))
for (g in unique(st$ko_gene)) {
  d <- st[st$ko_gene == g, ]
  fc <- d$FC
  n_sig <- sum(d$p.adj < 0.05)
  cat(sprintf("%-9s %7s %12.1f %12.3f %10d %10d\n",
              g, if (g %in% c("ACTA2", "YKT6", "TSHR")) "低" else "正常",
              max(fc), median(fc), sum(fc > 3.84), n_sig))
}

cat("\n===== 说明 =====\n")
cat("FC = distance^2 / mean(distance^2)；卡方检验下 FC>3.84 约对应单基因 p<0.05\n")
cat("若「FC>3.84 数」几乎只有 1（即被敲除基因自身），说明扰动没有传播到网络其他部分\n")

cat("\n===== 单个基因的距离分布（CSF1R）=====\n")
d <- st[st$ko_gene == "CSF1R", ]
v <- sort(d$distance, decreasing = TRUE)
cat("距离最大 5 个：", paste(round(v[1:5], 5), collapse = ", "), "\n")
cat("中位距离：", round(median(v), 6), "  最大/中位 =", round(v[1] / median(v), 1), "\n")
cat("非零距离个数：", sum(v > 0), "/", length(v), "\n")
cat("其中最大值对应的基因：", d$gene[which.max(d$distance)], "\n")
