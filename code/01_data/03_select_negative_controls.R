# 挑选能通过 scQC（表达率 >5%）的阴性对照基因
suppressPackageStartupMessages({ library(Matrix) })

X <- Matrix::readMM("work/data/counts.mtx")
genes <- readLines("work/data/genes.txt")
rownames(X) <- genes
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

rate <- Matrix::rowMeans(X != 0)
mean_expr <- Matrix::rowMeans(X)

## 候选：其他谱系标志基因 / 神经元、星形、少突胶质、内皮、外周免疫
cand <- c(
  "GFAP","AQP4","SLC1A2","SLC1A3",                 # 星形胶质
  "MBP","PLP1","MOG","OLIG1","OLIG2","CNP","MAG",  # 少突胶质
  "SNAP25","SYT1","RBFOX3","NEFL","NEFM","NEFH","DNM1", # 神经元
  "CLDN5","PECAM1","VWF","CDH5",                   # 内皮
  "COL1A1","COL1A2","DCN",                         # 成纤维
  "KRT8","KRT18",                                  # 上皮
  "ACTA2","MYH11","TAGLN",                         # 平滑肌
  "PTPRC","CD3E","CD19","NKG7","MS4A1","CD79A","IGHM", # 免疫谱系
  "HBB","HBA1","MYH7","EPCAM","TSHR","ALB"
)

cat("===== 阴性对照候选评估 =====\n")
cat(sprintf("%-9s %9s %11s %10s   %s\n", "基因", "表达率", "平均counts", "在矩阵", "评估"))
res <- data.frame()
for (g in cand) {
  if (g %in% genes) {
    r <- 100 * rate[g]
    verdict <- if (r >= 5 && r <= 45) "★ 可用（能进网络，表达中等）" else if (r < 5) "✗ 会被 scQC 过滤掉" else "△ 表达偏高"
    cat(sprintf("%-9s %8.2f%% %11.3f %10s   %s\n", g, r, mean_expr[g], "是", verdict))
    res <- rbind(res, data.frame(gene = g, rate = rate[g], mean = mean_expr[g]))
  } else {
    cat(sprintf("%-9s %9s %11s %10s   %s\n", g, "-", "-", "否", "不在矩阵中"))
  }
}

cat("\n===== 推荐阴性对照（按表达率 5%-45% 且平均 counts 最低排序）=====\n")
ok <- res[res$rate >= 0.05 & res$rate <= 0.45, ]
ok <- ok[order(ok$rate), ]
print(head(ok, 12), row.names = FALSE)
