# 检查批次是否与条件混淆
suppressPackageStartupMessages({ library(Matrix) })

meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)
X <- Matrix::readMM("work/data/counts.mtx")
lib <- Matrix::colSums(X)
names(lib) <- readLines("work/data/barcodes.txt")

cat("===== batch_sn x condition =====\n")
print(table(meta$batch_sn, meta$condition))

cat("\n===== batch_sn x lesion_type =====\n")
print(table(meta$batch_sn, meta$lesion_type))

cat("\n===== sample x batch =====\n")
print(table(meta$sample_id, meta$batch_sn))

cat("\n===== 每个批次的细胞数与文库量 =====\n")
tab <- do.call(rbind, lapply(split(meta, meta$batch_sn), function(d) {
  data.frame(batch = d$batch_sn[1], n_cells = nrow(d),
             conditions = paste(unique(d$condition), collapse = "/"),
             median_lib = round(median(lib[d$barcode])))
}))
print(tab, row.names = FALSE)

cat("\n===== CA vs CI 的样本深度对比（MS 内部比较）=====\n")
for (lt in c("CA", "CI")) {
  cells <- meta$barcode[meta$lesion_type == lt]
  q <- quantile(lib[cells], c(0, 0.25, 0.5, 0.75, 1))
  cat(sprintf("%-4s n=%5d  最小 %6.0f  中位 %6.0f  最大 %7.0f\n", lt, length(cells), q[1], q[3], q[5]))
}
