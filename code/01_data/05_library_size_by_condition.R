# 文库量阈值对两个条件的选择性影响
suppressPackageStartupMessages({ library(Matrix) })

X <- Matrix::readMM("work/data/counts.mtx")
rownames(X) <- readLines("work/data/genes.txt")
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)

lib <- Matrix::colSums(X)
names(lib) <- colnames(X)

cat("===== 每条件文库量分布 =====\n")
for (g in c("Control", "MS")) {
  cells <- meta$barcode[meta$condition == g]
  q <- quantile(lib[cells], c(0, 0.25, 0.5, 0.75, 1))
  cat(sprintf("%-8s n=%5d  最小 %6.0f  25%% %6.0f  中位 %6.0f  75%% %6.0f  最大 %7.0f\n",
              g, length(cells), q[1], q[2], q[3], q[4], q[5]))
}

cat("\n===== 不同文库量阈值下的保留细胞数 =====\n")
cat(sprintf("%10s %12s %12s %10s\n", "阈值", "Control", "MS", "MS/Control"))
for (th in c(200, 300, 500, 750, 1000, 1500)) {
  nc <- sum(lib[meta$barcode[meta$condition == "Control"]] > th)
  nm <- sum(lib[meta$barcode[meta$condition == "MS"]] > th)
  cat(sprintf("%10d %12d %12d %10.2f\n", th, nc, nm, nm / nc))
}

cat("\n===== 每个样本的细胞数与文库量中位数 =====\n")
tab <- do.call(rbind, lapply(split(meta, meta$sample_id), function(d) {
  data.frame(
    sample = d$sample_id[1],
    condition = d$condition[1],
    lesion = d$lesion_type[1],
    n_cells = nrow(d),
    median_lib = round(median(lib[d$barcode])),
    pass_1000 = sum(lib[d$barcode] > 1000)
  )
}))
print(tab[order(tab$condition, -tab$median_lib), ], row.names = FALSE)
