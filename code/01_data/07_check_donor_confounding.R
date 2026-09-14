# 核实病灶类型与患者是否重叠（判断 CA vs CI 是否被供体混淆）
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

cat("===== 患者 x 病灶类型（细胞数）=====\n")
print(table(meta$patient, meta$lesion_type))

cat("\n===== 患者 x 条件 =====\n")
print(table(meta$patient, meta$condition))

cat("\n===== 病灶类型各自的来源患者与占比 =====\n")
for (lt in c("CA", "CI", "Ctrl")) {
  sub <- meta[meta$lesion_type == lt, ]
  tb <- sort(table(sub$patient), decreasing = TRUE)
  cat(sprintf("\n[%s] 共 %d 细胞，来源 %d 位患者\n", lt, nrow(sub), length(tb)))
  print(round(100 * tb / sum(tb), 1))
}

ca_pat <- unique(meta$patient[meta$lesion_type == "CA"])
ci_pat <- unique(meta$patient[meta$lesion_type == "CI"])
ct_pat <- unique(meta$patient[meta$lesion_type == "Ctrl"])
cat("\n===== 患者集合重叠情况 =====\n")
cat("CA 患者:", paste(ca_pat, collapse = ", "), "\n")
cat("CI 患者:", paste(ci_pat, collapse = ", "), "\n")
cat("Ctrl 患者:", paste(ct_pat, collapse = ", "), "\n")
cat("CA ∩ CI:", paste(intersect(ca_pat, ci_pat), collapse = ", "), if (length(intersect(ca_pat, ci_pat)) == 0) " <== 完全不重叠" else "", "\n")
cat("MS ∩ Control:", paste(intersect(c(unique(meta$patient[meta$condition == "MS"])), unique(meta$patient[meta$condition == "Control"])), collapse = ", "), "\n")
