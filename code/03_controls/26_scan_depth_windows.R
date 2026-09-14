# 扫描不同深度窗口，寻找"细胞数足够 + 深度匹配良好"的平衡点
suppressPackageStartupMessages({ library(Matrix) })
X <- Matrix::readMM("work/data/counts.mtx")
colnames(X) <- readLines("work/data/barcodes.txt")
meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)
lib <- Matrix::colSums(X)
ct <- meta$barcode[meta$condition == "Control"]
ms <- meta$barcode[meta$condition == "MS"]

pick <- function(cells, lo, hi) {
  s <- cells[lib[cells] > lo & lib[cells] <= hi]
  l <- lib[s]
  s[!(l %in% grDevices::boxplot.stats(l)$out)]
}

cat(sprintf("%-16s %8s %8s %8s %8s %10s\n",
            "窗口", "Control", "MS", "中位CT", "中位MS", "KS p"))
for (w in list(c(500, Inf), c(800, 6000), c(1000, 5000), c(1500, 5000),
               c(1500, 4500), c(2000, 5000), c(2000, 4500))) {
  a <- pick(ct, w[1], w[2]); b <- pick(ms, w[1], w[2])
  if (length(a) < 100 || length(b) < 100) next
  ks <- suppressWarnings(stats::ks.test(lib[a], lib[b]))$p.value
  cat(sprintf("%-16s %8d %8d %8s %8s %10.2g\n",
              sprintf("%d-%s", w[1], ifelse(is.infinite(w[2]), "Inf", w[2])),
              length(a), length(b),
              format(round(median(lib[a])), big.mark = ","),
              format(round(median(lib[b])), big.mark = ","), ks))
}
