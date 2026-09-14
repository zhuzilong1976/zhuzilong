#!/usr/bin/env Rscript
# =============================================================================
# Stability checkpoints C3 and C4.
#
# 29_stability_analysis.R (03_controls) evaluates the first two pre-specified
# checkpoints and leaves C3 and C4 as pointers. This script completes them from
# the saved replicate objects:
#
#   C3  number of panel genes with FDR < 0.05 (empirical null) == 0 in every
#       replicate -- re-implements the test of 21_empirical_null_analysis.R
#       (Tables S4/S5) on each of the 20 replicate distance matrices.
#   C4  no gene with FDR < 0.05 in the paired condition comparison -- pairs
#       Control replicate i with MS replicate i and repeats the analysis of
#       24_condition_comparison.R on those two independently built networks.
#
# Inputs : work/stability/MSvC_{Control,MS}_rep01..10.rds   (28_stability_replicates.R)
# Outputs: outputs/results/vko_stability_checkpoints/
#            checkpoint_C3_empirical_null_by_replicate.csv
#            checkpoint_C4_paired_condition_comparison.csv
#
# Run from the repository root:  Rscript code/07_supplementary/67_stability_checkpoints_C3C4.R
# Consumed by 06_figures/54_figure_stability.R (Figure S2C, S2D).
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })

st_dir  <- "work/stability"
out_dir <- "outputs/results/vko_stability_checkpoints"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")
neg   <- c("ACTA2","TSHR","PECAM1","YKT6")

files <- sort(list.files(st_dir, pattern = "_rep[0-9]+\\.rds$", full.names = TRUE))
cat("replicate files:", length(files), "\n\n")

## ---------------------------------------------------------------- C3
cat("================ C3: panel genes under the empirical null ================\n")
c3 <- list()
for (f in files) {
  o <- readRDS(f)
  D <- o$distances
  od <- o$outdegree
  genes <- rownames(D)

  ## specificity z: standardise every target gene across all knockouts
  mu <- colMeans(D, na.rm = TRUE)
  s  <- apply(D, 2, stats::sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- NA
  Z <- sweep(sweep(D, 2, mu, "-"), 2, s, "/")
  Z[!is.finite(Z)] <- NA
  diag(Z) <- NA
  maxz <- apply(Z, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  maxz[!is.finite(maxz)] <- NA

  null_genes <- setdiff(genes, c(panel, neg))
  null_maxz  <- maxz[null_genes]
  tested     <- intersect(c(panel, neg), genes)
  p <- vapply(tested, function(g) {
    mz <- unname(maxz[g])
    if (is.na(mz)) return(NA_real_)
    (1 + sum(null_maxz >= mz, na.rm = TRUE)) / (1 + sum(!is.na(null_maxz)))
  }, numeric(1))
  padj <- stats::p.adjust(p, method = "fdr")
  cand <- panel[panel %in% tested]
  c3[[basename(f)]] <- data.frame(
    file = basename(f), tag = o$params$tag, rep = o$params$rep,
    n_genes = length(genes), n_null = sum(!is.na(null_maxz)),
    n_panel_tested = length(cand),
    n_candidate_fdr05 = sum(padj[names(padj) %in% cand] < 0.05, na.rm = TRUE),
    n_negative_fdr05  = sum(padj[names(padj) %in% neg]  < 0.05, na.rm = TRUE),
    min_padj_candidate = round(min(padj[names(padj) %in% cand], na.rm = TRUE), 3),
    top_candidate = names(padj)[which.min(ifelse(names(padj) %in% cand, padj, Inf))],
    stringsAsFactors = FALSE)
}
c3tab <- do.call(rbind, c3)
print(c3tab[, c("tag","rep","n_genes","n_null","n_candidate_fdr05",
                "n_negative_fdr05","min_padj_candidate","top_candidate")],
      row.names = FALSE)
write.csv(c3tab, file.path(out_dir, "checkpoint_C3_empirical_null_by_replicate.csv"),
          row.names = FALSE)
cat(sprintf("\nC3: replicates with >=1 candidate gene at FDR<0.05: %d / %d  -> %s\n",
            sum(c3tab$n_candidate_fdr05 > 0), nrow(c3tab),
            ifelse(all(c3tab$n_candidate_fdr05 == 0), "PASS", "FAIL")))
cat(sprintf("    smallest candidate FDR across replicates: %.3f\n",
            min(c3tab$min_padj_candidate)))

## ---------------------------------------------------------------- C4
cat("\n================ C4: paired condition comparison ================\n")
c4 <- list()
reps <- sort(unique(c3tab$rep))
for (r in reps) {
  fc <- file.path(st_dir, sprintf("MSvC_Control_rep%02d.rds", r))
  fm <- file.path(st_dir, sprintf("MSvC_MS_rep%02d.rds", r))
  if (!file.exists(fc) || !file.exists(fm)) next
  A <- as.matrix(readRDS(fc)$network)
  B <- as.matrix(readRDS(fm)$network)
  g <- intersect(rownames(A), rownames(B))
  A <- A[g, g]; B <- B[g, g]
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(A, B, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  DR <- DR[order(DR$p.value), ]
  c4[[as.character(r)]] <- data.frame(
    replicate = r, n_genes = nrow(DR),
    network_rho = round(stats::cor(as.vector(A), as.vector(B), method = "spearman"), 3),
    n_p05 = sum(DR$p.value < 0.05), n_fdr05 = sum(DR$p.adj < 0.05),
    min_padj = signif(min(DR$p.adj), 3), top_gene = DR$gene[1],
    stringsAsFactors = FALSE)
}
c4tab <- do.call(rbind, c4)
print(c4tab, row.names = FALSE)
write.csv(c4tab, file.path(out_dir, "checkpoint_C4_paired_condition_comparison.csv"),
          row.names = FALSE)
cat(sprintf("\nC4: replicate pairs with >=1 gene at FDR<0.05: %d / %d  -> %s\n",
            sum(c4tab$n_fdr05 > 0), nrow(c4tab),
            ifelse(all(c4tab$n_fdr05 == 0), "PASS", "FAIL")))
cat(sprintf("    smallest FDR across pairs: %.3g\n", min(c4tab$min_padj)))

cat("\nwritten to:", normalizePath(out_dir), "\n")
