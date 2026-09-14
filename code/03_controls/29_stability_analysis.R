#!/usr/bin/env Rscript
# =============================================================================
# Aggregate stability replicates and evaluate the pre-specified checkpoints.
#
# Reads  : work/stability/<tag>_rep*.rds  (produced by 28_stability_replicates.R)
# Writes : results/tables/Table_S10_stability_replicates.csv
#          results/tables/Table_S11_stability_summary.csv
#
# Checkpoints (pre-specified)
#   C1  containment rate per replicate >= 99% (allowing outdegree-0 genes)
#   C2  median Spearman rho of the per-gene max_z ranking between replicate
#       pairs >= 0.6
#   C3  the panel-vs-empirical-null conclusion (number of panel genes with
#       FDR < 0.05) is 0 in every replicate
#   C4  the direction of the condition comparison is unchanged: no gene reaches
#       FDR < 0.05 in any replicate
# =============================================================================
suppressPackageStartupMessages({ library(Matrix) })

st_dir <- "work/stability"
files <- list.files(st_dir, pattern = "_rep[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("no replicate files found in ", st_dir)

analysis <- function(o) {
  W <- o$network; D <- o$distances; od <- o$outdegree
  genes <- rownames(D)
  res <- do.call(rbind, lapply(genes, function(g) {
    d <- D[g, ]
    if (all(!is.finite(d))) return(NULL)
    FC <- d^2 / mean(d^2, na.rm = TRUE)
    p <- stats::pchisq(FC, df = 1, lower.tail = FALSE)
    padj <- stats::p.adjust(p, method = "fdr")
    sig <- names(padj)[padj < 0.05]
    targets <- colnames(W)[W[g, ] != 0]
    data.frame(ko = g, outdegree = od[g], n_sig = length(sig),
               all_inside = all(sig %in% union(targets, g)),
               mean_abs = mean(abs(d), na.rm = TRUE))
  }))
  ## specificity z across all knockouts of this replicate
  mu <- colMeans(D, na.rm = TRUE)
  s <- apply(D, 2, stats::sd, na.rm = TRUE); s[s == 0 | is.na(s)] <- NA
  Z <- sweep(sweep(D, 2, mu, "-"), 2, s, "/")
  Z[!is.finite(Z)] <- NA
  diag(Z) <- NA
  maxz <- apply(Z, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  list(res = res, maxz = maxz)
}

per_rep <- list()
for (f in files) {
  o <- readRDS(f)
  a <- analysis(o)
  r <- a$res
  neg <- r$outdegree == 0
  per_rep[[basename(f)]] <- data.frame(
    file = basename(f), tag = o$params$tag, rep = o$params$rep,
    n_genes = nrow(r), n_cells = o$params$n_cells,
    density = round(100 * mean(o$network != 0), 1),
    containment_all = round(100 * mean(r$all_inside), 2),
    containment_excl_deg0 = round(100 * mean(r$all_inside[!neg]), 2),
    n_deg0 = sum(neg), n_sig_total = sum(r$n_sig),
    max_z = round(max(a$maxz, na.rm = TRUE), 2),
    stringsAsFactors = FALSE)
  per_rep[[basename(f)]]$maxz_vec <- I(list(a$maxz))
}

tab <- do.call(rbind, per_rep)
vecs <- tab$maxz_vec
tab$maxz_vec <- NULL

## C2: pairwise Spearman rho of max_z vectors within the same tag
rhos <- c()
for (tg in unique(tab$tag)) {
  idx <- which(tab$tag == tg)
  if (length(idx) < 2) next
  for (i in seq_len(length(idx) - 1)) for (j in (i + 1):length(idx)) {
    a <- vecs[[idx[i]]]; b <- vecs[[idx[j]]]
    g <- intersect(names(a)[is.finite(a)], names(b)[is.finite(b)])
    if (length(g) < 10) next
    rhos <- c(rhos, suppressWarnings(stats::cor(a[g], b[g], method = "spearman")))
  }
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, "results/tables/Table_S10_stability_replicates.csv", row.names = FALSE)

summ <- data.frame(
  checkpoint = c("C1 containment >= 99% per replicate",
                 "C2 median pairwise rho >= 0.6",
                 "C3 n panel genes with FDR<0.05 == 0 in every replicate",
                 "C4 no gene with FDR<0.05 in the condition comparison"),
  value = c(sprintf("min %.2f%% (over %d replicates)", min(tab$containment_excl_deg0), nrow(tab)),
            if (length(rhos)) sprintf("median %.3f (range %.3f-%.3f, %d pairs)",
                                      median(rhos), min(rhos), max(rhos), length(rhos)) else "n/a",
            "see Table_S4/S5 per-replicate null analysis",
            "see per-replicate condition comparison"),
  pass = c(all(tab$containment_excl_deg0 >= 99),
           if (length(rhos)) median(rhos) >= 0.6 else NA,
           NA, NA))
utils::write.csv(summ, "results/tables/Table_S11_stability_summary.csv", row.names = FALSE)

cat("\n===== 稳定性汇总 =====\n")
print(tab[, c("tag", "rep", "n_cells", "density", "containment_excl_deg0",
              "n_deg0", "max_z")], row.names = FALSE)
if (length(rhos)) {
  cat(sprintf("\nC2：max_z 排名的重复间 Spearman rho 中位 %.3f（%d 对，范围 %.3f-%.3f）\n",
              median(rhos), length(rhos), min(rhos), max(rhos)))
}
cat("\n已写入 results/tables/Table_S10_stability_replicates.csv 与 Table_S11_stability_summary.csv\n")
