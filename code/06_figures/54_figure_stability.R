#!/usr/bin/env Rscript
# =============================================================================
# Figure S2: stability of the virtual-knockout outputs across independent cell
# subsamples (10 replicates per condition), i.e. the four pre-specified
# checkpoints C1-C4.
#
# Inputs : work/stability/*_rep*.rds              (produced by 28_…)
#          results/tables/Table_S10_stability_replicates.csv
#          outputs/results/vko_stability_checkpoints/checkpoint_C3_*.csv   (67_…)
#          outputs/results/vko_stability_checkpoints/checkpoint_C4_*.csv   (67_…)
# Output : results/figures/FigS2_stability.png (2400 px) and .pdf
#
# Note: the C3/C4 checkpoint files are produced by
# 07_supplementary/67_stability_checkpoints_C3C4.R, because
# 29_stability_analysis.R evaluates C1 and C2 only.
# =============================================================================
suppressPackageStartupMessages({ library(Matrix) })

st_dir  <- "work/stability"
ck_dir  <- "outputs/results/vko_stability_checkpoints"
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

col_ctl <- "#1f6fb4"
col_ms  <- "#c0392b"

## ---------------------------------------------------------------- inputs
tab <- utils::read.csv("results/tables/Table_S10_stability_replicates.csv",
                       stringsAsFactors = FALSE)

## pairwise Spearman rho of the max_z ranking, recomputed exactly as in
## 29_stability_analysis.R
maxz_of <- function(f) {
  D <- readRDS(f)$distances
  mu <- colMeans(D, na.rm = TRUE)
  s  <- apply(D, 2, stats::sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- NA
  Z <- sweep(sweep(D, 2, mu, "-"), 2, s, "/")
  Z[!is.finite(Z)] <- NA
  diag(Z) <- NA
  apply(Z, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
}
files <- sort(list.files(st_dir, pattern = "_rep[0-9]+\\.rds$", full.names = TRUE))
maxz  <- lapply(files, maxz_of)
names(maxz) <- sub("_rep[0-9]+\\.rds$", "", basename(files))
reps  <- as.integer(sub(".*_rep([0-9]+)\\.rds$", "\\1", basename(files)))

rho <- do.call(rbind, lapply(unique(names(maxz)), function(tg) {
  idx <- which(names(maxz) == tg)
  if (length(idx) < 2) return(NULL)
  out <- do.call(rbind, combn(idx, 2, simplify = FALSE, FUN = function(p) {
    a <- maxz[[p[1]]]; b <- maxz[[p[2]]]
    g <- intersect(names(a)[is.finite(a)], names(b)[is.finite(b)])
    data.frame(tag = tg, rho = suppressWarnings(
      stats::cor(a[g], b[g], method = "spearman")))
  }))
  out
}))

c3 <- try(utils::read.csv(file.path(ck_dir, "checkpoint_C3_empirical_null_by_replicate.csv"),
                          stringsAsFactors = FALSE), silent = TRUE)
c4 <- try(utils::read.csv(file.path(ck_dir, "checkpoint_C4_paired_condition_comparison.csv"),
                          stringsAsFactors = FALSE), silent = TRUE)
have34 <- !inherits(c3, "try-error") && !inherits(c4, "try-error")

lab <- ifelse(grepl("_Control_", tab$file), paste0("C", tab$rep), paste0("M", tab$rep))
col <- ifelse(grepl("_Control_", tab$file), col_ctl, col_ms)
o   <- order(grepl("_MS_", tab$file), tab$rep)

## ---------------------------------------------------------------- drawing
draw <- function() {
  op <- graphics::par(mfrow = c(2, 2), mar = c(4.6, 4.6, 3.2, 1.0),
                      oma = c(0, 0, 2.6, 0), family = "sans")

  ## A  C1 containment ------------------------------------------------------
  b <- graphics::barplot(tab$containment_excl_deg0[o], names.arg = lab[o],
                         col = col[o], border = NA, ylim = c(90, 106), las = 1,
                         cex.names = 0.72, ylab = "Containment (%)",
                         main = "A  C1 containment  \u2014  pass", cex.main = 1.15)
  graphics::abline(h = 99, lty = 2, col = "grey20")
  graphics::mtext("non-degenerate knockouts", side = 3, line = 0.1, cex = 0.68, col = "grey25")
  graphics::text(mean(b), 101.2,
                 labels = sprintf("100%% in %d / %d replicates\n(dashed line: threshold 99%%)",
                                  nrow(tab), nrow(tab)),
                 cex = 0.88, font = 2)

  ## B  C2 ranking correlation ---------------------------------------------
  h <- graphics::hist(rho$rho, breaks = seq(-0.5, 0.9, by = 0.1), plot = FALSE)
  graphics::plot(h, col = "grey80", border = "white", las = 1,
                 xlim = c(-0.55, 0.9), ylim = c(0, max(h$counts) * 1.25),
                 xlab = "", ylab = "replicate pairs",
                 main = "B  C2 ranking correlation  \u2014  fail",
                 cex.main = 1.15)
  graphics::mtext("Spearman rho between replicate max_z rankings",
                  side = 1, line = 2.6, cex = 0.8)
  graphics::abline(v = stats::median(rho$rho), lty = 1, lwd = 2, col = col_ctl)
  graphics::abline(v = 0.6, lty = 2, lwd = 2, col = "grey20")
  graphics::text(0.58, max(h$counts) * 1.2, "threshold 0.6", cex = 0.68, col = "grey20",
                 adj = c(1, 0))
  graphics::text(stats::median(rho$rho) - 0.03, max(h$counts) * 1.2,
                 sprintf("median %.3f", stats::median(rho$rho)),
                 cex = 0.72, col = col_ctl, font = 2, adj = c(1, 0))
  graphics::text(0.42, max(h$counts) * 0.62,
                 sprintf("%d pairs\n(from %d replicates)", nrow(rho), nrow(tab)),
                 cex = 0.78)

  if (have34) {
    ## C  C3 empirical null --------------------------------------------------
    c3o <- c3[order(grepl("_MS_", c3$file), c3$rep), ]
    graphics::barplot(c3o$min_padj_candidate,
                      names.arg = ifelse(grepl("_Control_", c3o$file),
                                         paste0("C", c3o$rep), paste0("M", c3o$rep)),
                      col = ifelse(grepl("_Control_", c3o$file), col_ctl, col_ms),
                      border = NA, ylim = c(0, 1.22), las = 1, cex.names = 0.72,
                      ylab = "smallest candidate FDR\nper replicate",
                      main = "C  C3 panel genes vs empirical null  \u2014  pass",
                      cex.main = 1.15)
    graphics::abline(h = 0.05, lty = 2, col = "grey20")
    graphics::text(mean(seq_along(c3o$rep) * 1.2 - 1.2), 1.04,
                   sprintf("no panel gene below FDR 0.05 in %d / %d replicates\n(dashed line: threshold 0.05)",
                           nrow(c3), nrow(c3)), cex = 0.85, font = 2)

    ## D  C4 paired comparison ----------------------------------------------
    graphics::barplot(c4$n_fdr05, names.arg = paste0("pair ", c4$replicate),
                      col = "grey45", border = NA, ylim = c(0, 3.6), las = 1,
                      cex.names = 0.72, ylab = "genes at FDR < 0.05",
                      yaxt = "n",
                      main = "D  C4 paired condition comparison  \u2014  fail",
                      cex.main = 1.15)
    graphics::axis(2, at = 0:3, las = 1)
    graphics::abline(h = 0, lty = 1, col = "grey20")
    for (i in seq_len(nrow(c4))) {
      graphics::text(i * 1.2 - 0.6, c4$n_fdr05[i] + 0.09,
                     ifelse(c4$n_fdr05[i] > 0, c4$top_gene[i], ""),
                     cex = 0.6, srt = 45)
    }
    graphics::text(nrow(c4) * 1.2 * 0.5, 3.3,
                   sprintf("primary analysis: 0 genes\nreplicates: %d / %d pairs with >=1",
                           sum(c4$n_fdr05 > 0), nrow(c4)), cex = 0.85, font = 2)
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "checkpoint_C3/C4 CSV not found\nsee code/07_supplementary/67_stability_checkpoints_C3C4.R",
                   cex = 0.9)
  }

  graphics::mtext(paste0("Stability of virtual-knockout outputs across independent cell subsamples ",
                         "(10 replicates per condition, 800 genes, 10 networks x 500 cells each)"),
                  side = 3, outer = TRUE, cex = 0.9, font = 2, family = "sans")
  graphics::par(op)
}

png(file.path(fig_dir, "FigS2_stability.png"), width = 2400, height = 1700, res = 190)
draw(); dev.off()
pdf(file.path(fig_dir, "FigS2_stability.pdf"), width = 12.6, height = 8.9)
draw(); dev.off()
cat("\nSaved results/figures/FigS2_stability.{png,pdf}\n")
cat(sprintf("C1: %s | C2 median rho %.3f | C3: %d sig in %d replicates | C4: %d/%d pairs with >=1 gene\n",
            ifelse(all(tab$containment_excl_deg0 >= 99), "pass", "fail"),
            stats::median(rho$rho),
            if (have34) sum(c3$n_candidate_fdr05 > 0) else NA, nrow(c3),
            if (have34) sum(c4$n_fdr05 > 0) else NA, nrow(c4)))
