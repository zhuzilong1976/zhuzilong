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
  op <- bib_par(mfrow = c(1, 3), mar = c(5.0, 4.8, 3.0, 0.8), oma = c(0, 0, 3.2, 0))

  ## A: affected gene counts
  vals <- sub$n_affected
  b <- graphics::barplot(vals, names.arg = sub$ko_gene, col = colors, border = NA,
                         ylab = "Genes with changed expression",
                         main = "A  KO affects many genes",
                         cex.names = BIB_CEX_AXIS)
  graphics::text(b, vals + max(vals) * 0.03, labels = vals, cex = BIB_CEX_TXT)
  graphics::mtext("n_propagation = 3", side = 3, line = 0.2, cex = BIB_CEX_TXT,
                  col = "grey30")

  ## B: locality - fraction of affected genes that are direct targets
  frac <- sub$pct_affected_direct
  lbl <- c(sub$ko_gene, "scTenifoldKnk\n(3,195 KOs)")
  val <- c(frac, 99.9)
  col <- c(colors, "grey25")
  b2 <- graphics::barplot(val, names.arg = lbl, col = col, border = NA, ylim = c(0, 115),
                          ylab = "Affected genes that are DIRECT targets (%)",
                          main = "B  Direct target share", cex.names = BIB_CEX_AXIS)
  graphics::abline(h = 100, lty = 2, col = "grey40")
  graphics::text(b2, val + 4, labels = sprintf("%.0f%%", val), cex = BIB_CEX_TXT)

  ## C: run-to-run variability of CellOracle
  run1 <- c(28, 1053, 732, 2230)     # first run
  run2 <- c(1205, 887, 2498, 2034)   # second run
  m <- rbind(run1, run2)
  colnames(m) <- sub$ko_gene
  b3 <- graphics::barplot(m, beside = TRUE, col = c("#95a5a6", "#34495e"), border = NA,
                          ylab = "Genes with changed expression",
                          main = "C  Run-to-run variation",
                          cex.names = BIB_CEX_AXIS,
                          legend.text = c("run 1", "run 2"),
                          args.legend = list(x = "topleft", bty = "n", cex = BIB_CEX_LEG))
  ## Kept inside the artwork because the 43-fold figure is not repeated in the
  ## manuscript text; the manuscript legend carries the rest of the panel text.
  graphics::mtext("Two independent CellOracle runs on the same data differ up to 43-fold",
                  side = 3, outer = TRUE, line = 2, cex = BIB_CEX_TXT, font = 2,
                  family = BIB_FONT)
  graphics::mtext("(SPI1: 28 vs 1,205); in scTenifoldKnk the affected set is deterministic given the network",
                  side = 3, outer = TRUE, line = 1, cex = BIB_CEX_TXT, font = 2,
                  family = BIB_FONT)
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_tool_comparison"),
           width_mm = BIB_FULL_MM, height_mm = 100)
cat("Saved results/figures/Fig_tool_comparison.{pdf,tiff,png}\n")
