# Supplementary figure: parameter sensitivity of the containment rule
suppressPackageStartupMessages({ library(Matrix) })
source("code/06_figures/_bib_figstyle.R")
psdir <- "outputs/results/param_sensitivity"
fig_dir <- "results/figures"

tags <- c("lam0_q090", "lam05_q090", "lam1_q090", "lam0_q095")
labels <- c("lambda=0\nq=0.90\n(default)", "lambda=0.5\nq=0.90",
            "lambda=1\nq=0.90", "lambda=0\nq=0.95")
## Two-line labels for the bars (the "(default)" wording stays in Table S2 and
## in the first setting's axis label, marked with an asterisk).
axis_labels <- c("lambda=0*\nq=0.90", "lambda=0.5\nq=0.90",
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
  ## Outer top margin 1.6 lines for the footnote that belongs to the whole
  ## figure; it used to sit in panel A's margin, where it read as a note on that
  ## panel alone. The panels keep the 8-line bottom margin their labels need.
  op <- bib_par(mfrow = c(1, 3), mar = c(8.0, 4.8, 3.0, 0.8), oma = c(0, 0, 1.6, 0))
  ## The manuscript writes these settings as nc_lambda / nc_q, so the axis uses
  ## the symbol rather than the word "lambda".
  ## lapply keeps these as language objects; as.expression turns the list into
  ## the expression vector barplot() needs to plot them as math.
  ax_lab <- as.expression(lapply(strsplit(axis_labels, "\n"), function(l) {
    bquote(lambda == .(sub("lambda=", "", l[1])) * "," ~ q == .(sub("q=", "", l[2])))
  }))

  ## A: compliance. All four bars sit near 100%, so the values go into the top
  ## margin where they cannot be clipped.
  b <- graphics::barplot(summ$compliance, names.arg = ax_lab, col = summ$col, border = NA,
                         ylim = c(0, 112), las = 2, cex.names = BIB_CEX_AXIS,
                         ylab = "Containment compliance (%)",
                         main = "A  Containment holds")
  graphics::abline(h = 100, lty = 2, col = "grey35")
  ## Rotated, like panels B and C: four horizontal labels at 8.1 pt are wider
  ## than the panel, so side-by-side they touch.
  graphics::text(b, summ$compliance + 3.5, sprintf("%.1f%%", summ$compliance),
                 cex = BIB_CEX_TXT, srt = 90, adj = c(0, 0.5))

  ## B: network density & median outdegree
  graphics::barplot(summ$density, names.arg = ax_lab, col = summ$col, border = NA,
                    ylim = c(0, 118), las = 2, cex.names = BIB_CEX_AXIS,
                    ylab = "Network density (%)",
                    main = "B  Density changes")
  graphics::text(b, summ$density + 4,
                 sprintf("%.1f%%\n(od=%d)", summ$density, round(summ$median_od)),
                 cex = BIB_CEX_TXT, srt = 90)

  ## C: number of outdegree-0 genes = the violations
  graphics::barplot(summ$n_deg0, names.arg = ax_lab, col = summ$col, border = NA,
                    ylim = c(0, max(summ$n_deg0) * 1.35), las = 2, cex.names = BIB_CEX_AXIS,
                    ylab = "Genes with outdegree 0",
                    main = "C  Outdegree-0 genes")
  graphics::text(b, summ$n_deg0 + max(summ$n_deg0) * 0.06,
                 labels = summ$n_deg0, cex = BIB_CEX_TXT, srt = 90)

  graphics::mtext("* default setting", side = 1, outer = TRUE, line = 0.4,
                  cex = BIB_CEX_TXT, col = "grey25")
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "FigS1_parameter_sensitivity"),
           width_mm = BIB_FULL_MM, height_mm = 108)
cat("\nSaved results/figures/FigS1_parameter_sensitivity.{pdf,tiff,png}\n")
