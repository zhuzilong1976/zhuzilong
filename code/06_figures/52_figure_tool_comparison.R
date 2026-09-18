# Figure: compare the locality of in-silico KO between scTenifoldKnk and CellOracle
source("code/06_figures/_bib_figstyle.R")
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
## The analysis output lives under outputs/ when CellOracle is re-run; the
## archived copy in results/tables/ is used when it is absent.
cmp_path <- c("outputs/results/celloracle/celloracle_vs_direct_targets.csv",
              "results/tables/Table_S9_celloracle_vs_direct_targets.csv")
cmp_path <- cmp_path[file.exists(cmp_path)][1]
cmp <- utils::read.csv(cmp_path, stringsAsFactors = FALSE)
sub <- cmp[cmp$threshold == "nonzero", ]

colors <- c("#c0392b", "#e67e22", "#2980b9", "#16a085")

draw <- function() {
  ## No outer banner: the two-line take-home message that used to run across the
  ## top of the artwork is now a one-line annotation inside panel C, which is
  ## where the run-to-run numbers are. mar top 3.6 keeps panel A's
  ## "n_propagation = 3" note clear of the plot box.
  ## mar right 2.2: panel C's title ("C  Run-to-run variation") is centred on the
  ## panel and its last characters reached the right edge of the device at 0.8.
  op <- bib_par(mfrow = c(1, 3), mar = c(5.0, 4.8, 3.6, 2.2), oma = c(0, 0, 0, 0))

  ## A: affected gene counts. The y-axis title matches panels B and C, and the
  ## panel title stays short: at 10.35 pt the longer "A  KO affects many genes"
  ## ran past the right edge of the panel and was cut off.
  vals <- sub$n_affected
  b <- graphics::barplot(vals, names.arg = sub$ko_gene, col = colors, border = NA,
                         ylim = c(0, max(vals) * 1.08),
                         ylab = "Genes with changed expression",
                         main = "A  KO-affected genes",
                         cex.names = BIB_CEX_AXIS)
  graphics::text(b, vals + max(vals) * 0.025, labels = vals, cex = BIB_CEX_TXT)
  graphics::mtext("n_propagation = 3", side = 3, line = 0.2, cex = BIB_CEX_TXT,
                  col = "grey30")

  ## B: locality - fraction of affected genes that are direct targets
  frac <- sub$pct_affected_direct
  lbl <- c(sub$ko_gene, "scTenifoldKnk\n(3,195 KOs)")
  val <- c(frac, 99.9)
  col <- c(colors, "grey25")
  ## Panel A already carries "Genes with changed expression" on its y axis.
  ## Repeating it here pushed the (long) title outside the panel, so the axis is
  ## labelled "DIRECT targets (%)" and the legend explains the rest.
  b2 <- graphics::barplot(val, names.arg = lbl, col = col, border = NA,
                          ylim = c(0, 112),
                          ylab = "DIRECT targets (%)",
                          main = "B  Direct target share", cex.names = BIB_CEX_AXIS)
  graphics::abline(h = 100, lty = 2, col = "grey40")
  graphics::text(b2, val + 4, labels = sprintf("%.0f%%", val), cex = BIB_CEX_TXT)

  ## C: run-to-run variability of CellOracle
  run1 <- c(28, 1053, 732, 2230)     # first run
  run2 <- c(1205, 887, 2498, 2034)   # second run
  m <- rbind(run1, run2)
  colnames(m) <- sub$ko_gene
  ## ylim has to be explicit. Left to itself barplot() sizes the axis for the
  ## first row only (max 2,230 here), so the 2,498 bar of run 2 was drawn
  ## outside the plot region and cut off at the panel edge. max(m) spans both
  ## rows and the 1.08 margin keeps the tallest bar clear of the top edge.
  b3 <- graphics::barplot(m, beside = TRUE, col = c("#95a5a6", "#34495e"), border = NA,
                          ylim = c(0, max(m) * 1.06),
                          ylab = "Genes with changed expression",
                          main = "C  Run-to-run variation",
                          cex.names = BIB_CEX_AXIS,
                          legend.text = c("run 1", "run 2"),
                          args.legend = list(x = "topleft", bty = "n", cex = BIB_CEX_LEG))
  ## Values on the bars: the 43-fold SPI1 contrast (28 vs 1,205) is the point of
  ## this panel, so it should not have to be read off the axis.
  graphics::text(b3, m + max(m) * 0.015, labels = m, cex = BIB_CEX_TXT)
  ## The 43-fold contrast does not appear in the manuscript text, so it stays in
  ## the artwork, as a single line above panel C. The range is 28 to 1,205 and 43
  ## is exactly round(1205 / 28).
  ## One short line, 25 characters: a longer version overran the 1.98 in panel on
  ## both sides, and a two-line version reached up into the panel title.
  graphics::mtext("SPI1 28 vs 1,205 (43-fold)",
                  side = 3, line = 0.3, cex = BIB_CEX_TXT, col = "grey25")
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_tool_comparison"),
           width_mm = BIB_FULL_MM, height_mm = 90)
cat("Saved results/figures/Fig_tool_comparison.{pdf,tiff,png}\n")
