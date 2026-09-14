# Figure: compare the locality of in-silico KO between scTenifoldKnk and CellOracle
fig_dir <- "outputs/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cmp <- utils::read.csv("outputs/results/celloracle/celloracle_vs_direct_targets.csv",
                       stringsAsFactors = FALSE)
sub <- cmp[cmp$threshold == "nonzero", ]

colors <- c("#c0392b", "#e67e22", "#2980b9", "#16a085")

draw <- function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(5.5, 4.6, 3.4, 1), oma = c(0, 0, 2.4, 0),
                      family = "sans")

  ## A: affected gene counts
  vals <- sub$n_affected
  b <- graphics::barplot(vals, names.arg = sub$ko_gene, col = colors, border = NA,
                         ylab = "Genes with changed expression",
                         main = "A  CellOracle: KO affects many genes", cex.main = 1.05,
                         cex.names = 0.85)
  graphics::text(b, vals + max(vals) * 0.03, labels = vals, cex = 0.8)
  graphics::mtext("n_propagation = 3", side = 3, line = -1.2, cex = 0.7, col = "grey30")

  ## B: locality - fraction of affected genes that are direct targets
  frac <- sub$pct_affected_direct
  lbl <- c(sub$ko_gene, "scTenifoldKnk\n(3,195 KOs)")
  val <- c(frac, 99.9)
  col <- c(colors, "grey25")
  b2 <- graphics::barplot(val, names.arg = lbl, col = col, border = NA, ylim = c(0, 115),
                          ylab = "Affected genes that are DIRECT targets (%)",
                          main = "B  One-hop vs multi-hop", cex.main = 1.05, cex.names = 0.8)
  graphics::abline(h = 100, lty = 2, col = "grey40")
  graphics::text(b2, val + 4, labels = sprintf("%.0f%%", val), cex = 0.8)

  ## C: run-to-run variability of CellOracle
  run1 <- c(28, 1053, 732, 2230)     # first run
  run2 <- c(1205, 887, 2498, 2034)   # second run
  m <- rbind(run1, run2)
  colnames(m) <- sub$ko_gene
  b3 <- graphics::barplot(m, beside = TRUE, col = c("#95a5a6", "#34495e"), border = NA,
                          ylab = "Genes with changed expression",
                          main = "C  CellOracle run-to-run variability",
                          cex.main = 1.05, cex.names = 0.85,
                          legend.text = c("run 1", "run 2"),
                          args.legend = list(x = "topleft", bty = "n", cex = 0.75))
  graphics::mtext(paste("Two independent CellOracle runs on the same data differ up to 43-fold (SPI1:",
                        "28 vs 1,205); in scTenifoldKnk the affected set is deterministic given the network"),
                  side = 3, outer = TRUE, cex = 0.72, font = 2, family = "sans")
  graphics::par(op)
}

png(file.path(fig_dir, "Fig_tool_comparison.png"), width = 2400, height = 950, res = 200)
draw(); dev.off()
pdf(file.path(fig_dir, "Fig_tool_comparison.pdf"), width = 12, height = 4.75)
draw(); dev.off()
cat("Saved outputs/figures/Fig_tool_comparison.{png,pdf}\n")
