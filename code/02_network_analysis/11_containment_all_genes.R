# 把"DR ⊆ {KO} ∪ 直接靶标"的验证扩展到全部基因（用已保存的 800x800 敲除距离矩阵）
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"

verify <- function(tag) {
  D <- readRDS(file.path(dir, paste0(tag, "_fullKO_distances.rds")))
  W <- as.matrix(readRDS(file.path(dir, paste0(tag, "_raw_result.rds")))$tensorNetworks$WT)
  od <- Matrix::rowSums(W != 0)
  genes <- rownames(D)

  res <- do.call(rbind, lapply(genes, function(g) {
    d <- D[g, ]
    if (all(!is.finite(d))) return(NULL)
    FC <- d^2 / mean(d^2, na.rm = TRUE)
    p <- stats::pchisq(FC, df = 1, lower.tail = FALSE)
    padj <- stats::p.adjust(p, method = "fdr")
    sig <- names(padj)[padj < 0.05]
    targets <- colnames(W)[W[g, ] != 0]
    allowed <- union(targets, g)
    data.frame(
      ko = g, outdegree = od[g], n_sig = length(sig),
      n_inside = sum(sig %in% allowed), all_inside = all(sig %in% allowed),
      frac_targets_sig = if (length(targets) > 0) length(intersect(sig, targets)) / length(targets) else NA,
      bound_ok = length(sig) <= od[g] + 1
    )
  }))
  res
}

for (tag in c("MSvC_Control", "MSvC_MS")) {
  r <- verify(tag)
  cat(sprintf("\n===== %s（n=%d 个基因）=====\n", tag, nrow(r)))
  cat(sprintf("显著集 ⊆ 靶标∪自身： %d / %d （%.1f%%）\n",
              sum(r$all_inside), nrow(r), 100 * mean(r$all_inside)))
  cat(sprintf("|DR| ≤ 出度+1：      %d / %d （%.1f%%）\n",
              sum(r$bound_ok), nrow(r), 100 * mean(r$bound_ok)))
  bad <- r[!r$all_inside, ]
  if (nrow(bad) > 0) {
    cat("违例基因（全部应为出度 0 的退化情形）：\n")
    print(head(bad[, c("ko", "outdegree", "n_sig", "n_inside")], 10), row.names = FALSE)
    cat(sprintf("其中出度为 0 的占：%d / %d\n", sum(bad$outdegree == 0), nrow(bad)))
  }
  cat(sprintf("\n出度分层：\n"))
  r$bin <- cut(r$outdegree, breaks = c(-1, 0, 100, 250, 400, 550, 700, 1000))
  print(do.call(rbind, lapply(split(r, r$bin), function(x) data.frame(
    outdegree_bin = as.character(x$bin[1]), n = nrow(x),
    median_n_sig = median(x$n_sig),
    median_frac_targets_sig = round(median(x$frac_targets_sig, na.rm = TRUE), 3),
    pct_all_inside = round(100 * mean(x$all_inside), 1)
  ))), row.names = FALSE)
  utils::write.csv(r, file.path(dir, paste0("containment_", sub("MSvC_", "", tag), ".csv")),
                   row.names = FALSE)
}
