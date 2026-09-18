#!/usr/bin/env Rscript
# =============================================================================
# Mechanistic test of weight dilution (answers reviewer R2).
#
# For a knockout g and one of its direct targets j, define the incoming-weight
# share  s(g,j) = W[g,j] / sum_k W[k,j]  (rows of W are regulators). The dilution
# account predicts that j is called differentially regulated in proportion to
# s(g,j), and that no target reaches significance once the knockout's share is
# diluted by hundreds of other regulators.
#
# Knockouts are sampled across the outdegree range of each network (three near
# each of five quantiles, outdegree-0 genes excluded because their knockout is a
# no-op). Significance uses the empirical null, as in the panel analysis
# (10_run_vko_panel.R with empirical_null = TRUE), so that enough targets are
# positive to test the relationship.
#
# Inputs : outputs/results/vko_ms_microglia/MSvC_{Control,MS}_raw_result.rds
#          outputs/results/vko_ms_microglia/MSvC_{Control,MS}_fullKO_distances.rds
# Outputs: results/tables/Table_S14_weight_share_vs_significance.csv
#          results/figures/FigS3_weight_share.{pdf,png,tiff}
# Run from the repository root.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(scTenifoldKnk)
})
source("code/06_figures/_bib_figstyle.R")

dir  <- "outputs/results/vko_ms_microglia"
tags <- c(Control = "MSvC_Control", MS = "MSvC_MS")

rows <- list()
for (nm in names(tags)) {
  tag <- tags[[nm]]
  W <- as.matrix(readRDS(file.path(dir, paste0(tag, "_raw_result.rds")))$tensorNetworks$WT)
  D <- as.matrix(readRDS(file.path(dir, paste0(tag, "_fullKO_distances.rds"))))
  g <- intersect(rownames(W), rownames(D))
  W <- W[g, g, drop = FALSE]; D <- D[g, g, drop = FALSE]
  od  <- Matrix::rowSums(W != 0)
  inc <- Matrix::colSums(W)
  inc[inc <= 0] <- NA_real_

  ## Outdegree is strongly skewed in these networks (most genes exceed 500
  ## targets), so quantile sampling would land entirely in the diluted region.
  ## Sample the low tail explicitly, plus the median and the upper decile.
  cand <- g[od > 0]
  ord  <- order(od[cand])
  low  <- cand[ord][seq_len(min(12, length(cand)))]
  med  <- cand[order(abs(od[cand] - stats::median(od[cand])))][1:4]
  hi   <- cand[order(abs(od[cand] - stats::quantile(od[cand], 0.95)))][1:4]
  sel  <- unique(c(low, med, hi))
  cat(sprintf("%s: %d knockouts sampled, outdegree %.0f-%.0f (median %.0f)\n",
              nm, length(sel), min(od[sel]), max(od[sel]), stats::median(od[cand])))

  for (ko in sel) {
    KO <- W; KO[ko, ] <- 0
    MA <- scTenifoldNet::manifoldAlignment(W, KO, d = 2)
    DR <- scTenifoldKnk::dRegulation(MA, empiricalNull = TRUE)
    sig_genes <- DR$gene[DR$p.adj < 0.05]

    tgt <- setdiff(g[W[ko, ] != 0], ko)
    if (length(tgt) == 0) next
    rows[[paste(tag, ko)]] <- data.frame(
      network = nm, knockout = ko, outdegree = as.integer(od[ko]),
      target = tgt, share = as.numeric(W[ko, tgt] / inc[tgt]),
      significant = tgt %in% sig_genes, stringsAsFactors = FALSE)
  }
}
D <- do.call(rbind, rows)
ro <- D[is.finite(D$share), ]
ro$bin <- cut(ro$share, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)

## pooled relationship -------------------------------------------------------
binned <- do.call(rbind, lapply(split(ro, ro$bin), function(x) {
  data.frame(bin = as.character(x$bin[1]), share_mid = mean(x$share),
             n = nrow(x), pct_significant = round(100 * mean(x$significant), 1))
}))
binned <- binned[order(binned$share_mid), ]

## effect size: rank-based AUC and per-knockout Spearman ---------------------
auc <- function(x, y) {
  r <- rank(x); n1 <- sum(y); n0 <- sum(!y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
overall_auc <- auc(ro$share, ro$significant)
per_ko <- do.call(rbind, lapply(split(ro, paste(ro$network, ro$knockout)), function(x) {
  rho <- if (length(unique(x$significant)) > 1)
    suppressWarnings(stats::cor(x$share, as.numeric(x$significant), method = "spearman"))
  else NA_real_
  data.frame(network = x$network[1], knockout = x$knockout[1],
             outdegree = x$outdegree[1], n_targets = nrow(x),
             n_significant = sum(x$significant),
             median_share_sig = round(stats::median(x$share[x$significant]), 3),
             median_share_non = round(stats::median(x$share[!x$significant]), 3),
             rho = round(rho, 3))
}))
per_ko <- per_ko[order(per_ko$network, per_ko$outdegree), ]

utils::write.csv(binned, "results/tables/Table_S14_weight_share_vs_significance.csv",
                 row.names = FALSE)
utils::write.csv(per_ko, "results/tables/Table_S14b_weight_share_per_knockout.csv",
                 row.names = FALSE)
utils::write.csv(ro, "results/tables/Table_S14c_weight_share_target_level.csv",
                 row.names = FALSE)

cat("\n===== pooled share vs significance =====\n")
print(binned, row.names = FALSE)
cat(sprintf("\nrank AUC (share predicts significance): %.3f | target pairs: %d | positives: %d\n",
            overall_auc, nrow(ro), sum(ro$significant)))
cat("\n===== per knockout =====\n")
print(per_ko, row.names = FALSE)

## figure --------------------------------------------------------------------
## Panel B uses the knockouts whose targets are all significant, to show that the
## per-target share does not decide the call within such a knockout. (Within that
## subset every target is positive, so a rank AUC is undefined there; the pooled
## AUC across all sampled knockouts is reported instead.)
low <- ro[ro$network == "Control" & ro$knockout %in% per_ko$knockout[per_ko$n_significant > 0], ]
share_rng <- range(low$share)

draw <- function() {
  ## mar right 2.0: panel B's x-axis title is longer than the panel itself, so at
  ## 0.8 the last characters were cut by the right edge of the device.
  op <- bib_par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.4, 2.4))

  cols <- ifelse(per_ko$network == "Control", BIB_TINT$blue, BIB_TINT$red)
  graphics::plot(per_ko$outdegree, 100 * per_ko$n_significant / per_ko$n_targets,
                 log = "x", pch = 21, bg = cols, col = "black", cex = 1.4,
                 xlab = "outdegree of the knocked-out gene (log scale)",
                 ylab = "direct targets called significant (%)",
                 main = "A  All-or-none switch with outdegree", las = 1,
                 ylim = c(-10, 112), xlim = c(0.8, 2000))
  graphics::abline(h = c(0, 100), lty = 2, col = "grey45")
  graphics::legend("bottomleft", legend = c("Control network", "MS network"),
                   pch = 21, pt.bg = c(BIB_TINT$blue, BIB_TINT$red),
                   col = "black", bty = "n", cex = BIB_CEX_LEG, pt.cex = 1.2)

  h <- graphics::hist(low$share, breaks = 30, plot = FALSE)
  graphics::plot(h, col = BIB_TINT$blue, border = "white", las = 1,
                 xlim = c(0, max(0.07, max(low$share))),
                 xlab = "incoming-weight share from the knockout",
                 ylab = "targets called significant",
                 main = "B  Significant targets carry minute shares")
  ## One short line: the earlier single line was wider than the panel and was cut
  ## by the device edge, and a stacked version reached into the panel title. The
  ## outdegree range and the largest share are given in the figure legend.
  graphics::mtext(sprintf("all %d targets significant", nrow(low)),
                  side = 3, line = 0.25, cex = BIB_CEX_TXT, col = "grey25")
  graphics::abline(v = median(low$share), lty = 2, col = "grey35")
  graphics::text(median(low$share), max(h$counts) * 0.9,
                 sprintf("median %.3f", median(low$share)),
                 pos = 4, cex = BIB_CEX_TXT, col = "grey25")
  graphics::par(op)
}

bib_render(draw, "results/figures/FigS3_weight_share", width_mm = BIB_FULL_MM, height_mm = 82)
cat("\nwritten: results/tables/Table_S14*.csv, results/figures/FigS3_weight_share.{pdf,png,tiff}\n")
