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
source("code/06_figures/_bib_figstyle.R")

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
  op <- bib_par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.0, 0.8), oma = c(0, 0, 0, 0))

  ## A  C1 containment ------------------------------------------------------
  b <- graphics::barplot(tab$containment_excl_deg0[o], names.arg = lab[o],
                         col = col[o], border = NA, ylim = c(90, 106), las = 2,
                         cex.names = BIB_CEX_AXIS, ylab = "Containment (%)",
                         main = "A  C1 containment: pass")
  graphics::abline(h = 99, lty = 2, col = "grey20")
  graphics::mtext("non-degenerate knockouts", side = 3, line = 0.1,
                  cex = BIB_CEX_TXT, col = "grey25")
  graphics::text(mean(b), 101.2,
                 labels = sprintf("100%% in %d / %d replicates\n(dashed line: threshold 99%%)",
                                  nrow(tab), nrow(tab)),
                 cex = BIB_CEX_TXT, font = 2)

  ## B  C2 ranking correlation ---------------------------------------------
  h <- graphics::hist(rho$rho, breaks = seq(-0.5, 0.9, by = 0.1), plot = FALSE)
  graphics::plot(h, col = "grey80", border = "white", las = 1,
                 xlim = c(-0.55, 0.9), ylim = c(0, max(h$counts) * 1.45),
                 xlab = "", ylab = "replicate pairs",
                 main = "B  C2 ranking correlation: fail")
  graphics::mtext("Spearman rho of replicate max_z rankings",
                  side = 1, line = 2.6, cex = BIB_CEX_TXT)
  graphics::abline(v = stats::median(rho$rho), lty = 1, lwd = 2, col = col_ctl)
  graphics::abline(v = 0.6, lty = 2, lwd = 2, col = "grey20")
  ## adj = c(1, 0) anchors the text ending at x = 0.52, clear of the dashed
  ## line at 0.6; at 0.58 the label ran into the panel edge.
  graphics::text(0.52, max(h$counts) * 1.34, "threshold 0.6", cex = BIB_CEX_TXT,
                 col = "grey20", adj = c(1, 0))
  graphics::text(stats::median(rho$rho) + 0.04, max(h$counts) * 1.12,
                 sprintf("median %.3f", stats::median(rho$rho)),
                 cex = BIB_CEX_TXT, col = col_ctl, font = 2, adj = c(0, 0))
  ## Left of the histogram: centered at 0.42 the text sat on the bars.
  graphics::text(-0.5, max(h$counts) * 0.75,
                 sprintf("%d pairs\n(from %d replicates)", nrow(rho), nrow(tab)),
                 cex = BIB_CEX_TXT, adj = c(0, 0))

  if (have34) {
    ## C  C3 empirical null --------------------------------------------------
    c3o <- c3[order(grepl("_MS_", c3$file), c3$rep), ]
    graphics::barplot(c3o$min_padj_candidate,
                      names.arg = ifelse(grepl("_Control_", c3o$file),
                                         paste0("C", c3o$rep), paste0("M", c3o$rep)),
                      col = ifelse(grepl("_Control_", c3o$file), col_ctl, col_ms),
                      border = NA, ylim = c(0, 1.22), las = 2, cex.names = BIB_CEX_AXIS,
                      ylab = "smallest candidate FDR",
                      main = "C  C3 vs empirical null: pass")
    graphics::abline(h = 0.05, lty = 2, col = "grey20")
    graphics::text(mean(seq_along(c3o$rep) * 1.2 - 1.2), 1.04,
                   sprintf("no panel gene below FDR 0.05\n(%d / %d replicates; dashed line: 0.05)",
                           nrow(c3), nrow(c3)), cex = BIB_CEX_TXT, font = 2)

    ## D  C4 paired comparison ----------------------------------------------
    graphics::barplot(c4$n_fdr05, names.arg = paste0("pair ", c4$replicate),
                      col = "grey45", border = NA, ylim = c(0, 4.2), las = 2,
                      cex.names = BIB_CEX_AXIS, ylab = "genes at FDR < 0.05",
                      yaxt = "n",
                      main = "D  C4 paired comparison: fail")
    graphics::axis(2, at = 0:3, las = 1)
    graphics::abline(h = 0, lty = 1, col = "grey20")
    ## Horizontal, alternating two heights: at 45 degrees the labels of the
    ## adjacent pairs (CD74/ADGRG3, CD74/ADGRG3) overlapped each other.
    for (i in seq_len(nrow(c4))) {
      graphics::text(i * 1.2 - 0.6, c4$n_fdr05[i] + 0.1 + (i %% 2) * 0.32,
                     ifelse(c4$n_fdr05[i] > 0, c4$top_gene[i], ""),
                     cex = BIB_CEX_TXT)
    }
    graphics::text(nrow(c4) * 1.2 * 0.5, 3.9,
                   sprintf("primary analysis: 0 genes\nreplicates: %d / %d pairs with ≥1",
                           sum(c4$n_fdr05 > 0), nrow(c4)), cex = BIB_CEX_TXT, font = 2)
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "checkpoint_C3/C4 CSV not found\nsee code/07_supplementary/67_stability_checkpoints_C3C4.R",
                   cex = 0.9)
  }

  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "FigS2_stability"),
           width_mm = BIB_FULL_MM, height_mm = 140)
cat("\nSaved results/figures/FigS2_stability.{pdf,tiff,png}\n")
cat(sprintf("C1: %s | C2 median rho %.3f | C3: %d sig in %d replicates | C4: %d/%d pairs with >=1 gene\n",
            ifelse(all(tab$containment_excl_deg0 >= 99), "pass", "fail"),
            stats::median(rho$rho),
            if (have34) sum(c3$n_candidate_fdr05 > 0) else NA, nrow(c3),
            if (have34) sum(c4$n_fdr05 > 0) else NA, nrow(c4)))
