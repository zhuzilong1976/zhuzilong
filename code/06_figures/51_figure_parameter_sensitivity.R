# Supplementary figure: parameter sensitivity of the containment rule
suppressPackageStartupMessages({ library(Matrix) })
psdir <- "outputs/results/param_sensitivity"
fig_dir <- "outputs/figures"

tags <- c("lam0_q090", "lam05_q090", "lam1_q090", "lam0_q095")
labels <- c("lambda=0\nq=0.90\n(default)", "lambda=0.5\nq=0.90",
            "lambda=1\nq=0.90", "lambda=0\nq=0.95")
cols <- c("#1f6fb4", "#5dade2", "#7d3c98", "#c0392b")

summ <- do.call(rbind, lapply(seq_along(tags), function(i) {
  f <- file.path(psdir, paste0(tags[i], ".rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  W <- o$net; res <- o$res
  od <- Matrix::rowSums(W != 0)
  data.frame(tag = tags[i], label = labels[i], col = cols[i],
             density = 100 * mean(W != 0),
             median_od = median(od),
             n_deg0 = sum(od == 0),
             n = nrow(res),
             compliance = 100 * mean(res$all_inside),
             median_nsig = median(res$n_sig))
}))
utils::write.csv(summ[, setdiff(names(summ), "col")],
                 file.path(fig_dir, "param_sensitivity_summary.csv"), row.names = FALSE)
print(summ[, c("label", "density", "median_od", "n_deg0", "compliance", "median_nsig")],
      row.names = FALSE)

draw <- function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(6.5, 4.6, 3.2, 1), oma = c(0, 0, 2.4, 0),
                      family = "sans")

  ## A: compliance
  b <- graphics::barplot(summ$compliance, names.arg = summ$label, col = summ$col, border = NA,
                         ylim = c(0, 108), las = 1, cex.names = 0.7,
                         ylab = "Containment compliance (%)",
                         main = "A  Rule holds under all settings", cex.main = 1.05)
  graphics::abline(h = 100, lty = 2, col = "grey35")
  graphics::text(b, summ$compliance + 3,
                 labels = sprintf("%.1f%%", summ$compliance), cex = 0.8)

  ## B: network density & median outdegree
  graphics::barplot(summ$density, names.arg = summ$label, col = summ$col, border = NA,
                    ylim = c(0, 105), las = 1, cex.names = 0.7,
                    ylab = "Network density (%)",
                    main = "B  Parameters do change the network", cex.main = 1.05)
  graphics::text(b, summ$density + 3, sprintf("%.1f%%\n(od=%d)", summ$density, round(summ$median_od)),
                 cex = 0.75)

  ## C: number of outdegree-0 genes = the violations
  graphics::barplot(summ$n_deg0, names.arg = summ$label, col = summ$col, border = NA,
                    ylim = c(0, max(summ$n_deg0) * 1.35), las = 1, cex.names = 0.7,
                    ylab = "Genes with outdegree 0",
                    main = "C  All violations are outdegree-0 genes", cex.main = 1.05)
  graphics::text(b, summ$n_deg0 + max(summ$n_deg0) * 0.06,
                 labels = summ$n_deg0, cex = 0.9)
  graphics::mtext(paste("Parameter sensitivity of the containment rule (400-gene control network,",
                        "10 networks x 400 cells; n = 403 genes)"),
                  side = 3, outer = TRUE, cex = 0.82, font = 2, family = "sans")
  graphics::par(op)
}

png(file.path(fig_dir, "FigS1_parameter_sensitivity.png"), width = 2400, height = 950, res = 200)
draw(); dev.off()
pdf(file.path(fig_dir, "FigS1_parameter_sensitivity.pdf"), width = 12, height = 4.75)
draw(); dev.off()
cat("\nSaved outputs/figures/FigS1_parameter_sensitivity.{png,pdf}\n")
