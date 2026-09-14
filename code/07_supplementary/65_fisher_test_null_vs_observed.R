m <- matrix(c(58, 800 - 58, 45, 797 - 45), nrow = 2, byrow = TRUE)
rownames(m) <- c("observed_Control_vs_MS", "null_halfA_vs_halfB")
colnames(m) <- c("p_lt_0.05", "not_significant")
print(m)
ft <- stats::fisher.test(m)
cat(sprintf("\nFisher 精确检验 p = %.3f\n", ft$p.value))
cat(sprintf("比值比 = %.3f，95%% CI = %.3f – %.3f\n", ft$estimate, ft$conf.int[1], ft$conf.int[2]))
cat(sprintf("\n观测显著率 %.2f%% vs 对照假阳性率 %.2f%%，绝对差 %.2f 个百分点\n",
            100 * 58 / 800, 100 * 45 / 797, 100 * 58 / 800 - 100 * 45 / 797))
