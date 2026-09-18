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
  ## Two rows, not one: in a 1 x 3 layout each panel was about 40 mm wide, which
  ## is too narrow for four x-axis labels (SPI1, IRF8, STAT3, NFKB1 are 8 to 14 mm
  ## each) and far too narrow for the eight bars of panel C. A and B now share the
  ## top row at about 78 mm each, and C spans the full 157 mm underneath, so the
  ## bars are wide and every value label has room.
  ##
  ## No outer banner either: the two-line take-home message that used to run
  ## across the top is panel C's one-line annotation, where the run-to-run numbers
  ## are. mar top 3.6 keeps panel A's "n_propagation = 3" note clear of the box.
  op <- bib_par(mar = c(5.0, 4.8, 3.6, 1.6), oma = c(0, 0, 0, 0))
  graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))

  ## A: affected gene counts. The y-axis title matches panels B and C, and the
  ## panel title stays short: at 10.35 pt the longer "A  KO affects many genes"
  ## ran past the right edge of the panel and was cut off.
  vals <- sub$n_affected
  b <- graphics::barplot(vals, names.arg = sub$ko_gene, col = colors, border = NA,
                         ylim = c(0, max(vals) * 1.08),
                         ylab = "Genes with changed expression",
                         main = "A  KO-affected genes",
                         cex.names = BIB_CEX_AXIS)
  ## adj = c(0.5, 0): the text sits on top of the bar instead of straddling its
  ## top edge, which is what the default centred adj did.
  graphics::text(b, vals + max(vals) * 0.02, labels = vals, cex = BIB_CEX_TXT,
                 adj = c(0.5, 0))
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
  graphics::text(b2, val + 2.5, labels = sprintf("%.0f%%", val), cex = BIB_CEX_TXT,
                 adj = c(0.5, 0))

  ## C: run-to-run variability of CellOracle
  run1 <- c(28, 1053, 732, 2230)     # first run
  run2 <- c(1205, 887, 2498, 2034)   # second run
  m <- rbind(run1, run2)
  colnames(m) <- sub$ko_gene
  ## The value labels need headroom above the tallest bar, and ylim has to be
  ## explicit anyway: left to itself barplot() sizes the axis for the first row
  ## only (max 2,230), so the 2,498 bar of run 2 was drawn outside the plot
  ## region. max(m) spans both rows and 1.14 leaves room for the labels.
  graphics::par(mar = c(4.6, 4.8, 3.4, 1.6))
  b3 <- graphics::barplot(m, beside = TRUE, col = c("#95a5a6", "#34495e"), border = NA,
                          ylim = c(0, max(m) * 1.14),
                          ylab = "Genes with changed expression",
                          main = "C  Run-to-run variation",
                          cex.names = BIB_CEX_AXIS,
                          legend.text = c("run 1", "run 2"),
                          args.legend = list(x = "topleft", bty = "n", cex = BIB_CEX_LEG))
  ## Values on the bars: the 43-fold SPI1 contrast (28 vs 1,205) is the point of
  ## this panel, so it should not have to be read off the axis. adj = c(0.5, 0)
  ## puts each number above its bar rather than across the bar top, and the extra
  ## headroom in ylim keeps all eight of them inside the panel.
  graphics::text(b3, m + max(m) * 0.02, labels = m, cex = BIB_CEX_TXT,
                 adj = c(0.5, 0))
  ## The 43-fold contrast does not appear in the manuscript text, so it stays in
  ## the artwork, as a single line above panel C. The range is 28 to 1,205 and 43
  ## is exactly round(1205 / 28).
  graphics::mtext("SPI1 28 vs 1,205 (43-fold)",
                  side = 3, line = 0.3, cex = BIB_CEX_TXT, col = "grey25")
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_tool_comparison"),
           width_mm = BIB_FULL_MM, height_mm = 150)
cat("Saved results/figures/Fig_tool_comparison.{pdf,tiff,png}\n")
