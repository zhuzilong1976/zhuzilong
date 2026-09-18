# Publication-quality figure: the DR gene set is determined by the KO gene's outdegree
# All labels in English to avoid font-rendering issues on Windows
suppressPackageStartupMessages({ library(Matrix) })
source("code/06_figures/_bib_figstyle.R")
dir <- "outputs/results/vko_ms_microglia"
ndir <- "outputs/results/technical_null"
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

containment <- function(D, W) {
  od <- Matrix::rowSums(W != 0)
  genes <- rownames(D)
  do.call(rbind, lapply(genes, function(g) {
    d <- D[g, ]
    if (all(!is.finite(d))) return(NULL)
    FC <- d^2 / mean(d^2, na.rm = TRUE)
    p <- stats::pchisq(FC, df = 1, lower.tail = FALSE)
    padj <- stats::p.adjust(p, method = "fdr")
    sig <- names(padj)[padj < 0.05]
    targets <- colnames(W)[W[g, ] != 0]
    data.frame(ko = g, outdegree = od[g], n_sig = length(sig),
               all_inside = all(sig %in% union(targets, g)),
               frac_targets_sig = if (length(targets) > 0)
                 length(intersect(sig, targets)) / length(targets) else NA)
  }))
}

specs <- list(
  list(nm = "MS microglia, control net", D = file.path(dir, "MSvC_Control_fullKO_distances.rds"),
       Wr = file.path(dir, "MSvC_Control_raw_result.rds"), nested = TRUE, col = "#1f6fb4"),
  list(nm = "MS microglia, lesion net", D = file.path(dir, "MSvC_MS_fullKO_distances.rds"),
       Wr = file.path(dir, "MSvC_MS_raw_result.rds"), nested = TRUE, col = "#c0392b"),
  list(nm = "Control split-half A", D = file.path(ndir, "CtrlHalfA_fullKO.rds"),
       Wr = file.path(ndir, "CtrlHalfA_network.rds"), nested = FALSE, col = "#7fb3d5"),
  list(nm = "Control split-half B", D = file.path(ndir, "CtrlHalfB_fullKO.rds"),
       Wr = file.path(ndir, "CtrlHalfB_network.rds"), nested = FALSE, col = "#e59866")
)

res <- list()
for (s in specs) {
  if (!file.exists(s$D) || !file.exists(s$Wr)) { cat("skip:", s$nm, "\n"); next }
  W <- readRDS(s$Wr)
  if (s$nested) W <- W$tensorNetworks$WT
  W <- as.matrix(W)
  r <- containment(readRDS(s$D), W)
  r$network <- s$nm; r$col <- s$col
  res[[s$nm]] <- r
  cat(sprintf("%-30s n=%d  containment %.1f%%\n", s$nm, nrow(r), 100 * mean(r$all_inside)))
}
all_res <- do.call(rbind, res)
utils::write.csv(all_res, file.path(fig_dir, "containment_all_networks.csv"), row.names = FALSE)

## The developers' published Trem2 result. The RDS written by
## code/04_replication/33_author_trem2_validation.R is not kept in the repo
## (it is derived from GPL-3 author data), so the two reported values are used
## directly when it is absent.
aut_path <- "work/author/author_validation.rds"
if (file.exists(aut_path)) {
  aut <- readRDS(aut_path)
  aut_od <- as.numeric(aut$od[aut$gKO])
  aut_n <- aut$n_sig
} else {
  aut_od <- 2744
  aut_n  <- 128
  cat("author_validation.rds not found; using published values od=2744, |DR|=128\n")
}

## Identity colours, kept identical in every panel (points, bars, legend).
cols <- vapply(res, function(x) x$col[1], character(1))

draw <- function() {
  ## The x-axis title is one short string shared by panels A and B, so the
  ## column gap can stay small. A longer title would run past the left edge of
  ## the device, because it is centred on the panel and panel A is flush left.
  ## Bottom margin 6.0 and top 3.4: panel C carries two-line rotated network
  ## labels under the axis and a row of rotated value labels above the bars, and
  ## at mar = 4.6 / 3.0 both ran past the device edge and were cut (measured as
  ## ink on the last raster row of the exported PNG).
  op <- bib_par(mfrow = c(1, 3), mar = c(8.6, 4.4, 3.4, 0.8), oma = c(0, 0, 0, 0))
  xlab_ko <- "Outdegree of KO gene"

  ## Panel A
  graphics::plot(NA, xlim = c(1, 20000), ylim = c(0.8, 3000), log = "xy",
                 xlab = xlab_ko,
                 ylab = "absolute DR genes",
                 main = "A  |DR| vs outdegree")
  graphics::abline(0, 1, lty = 2, col = "grey35", lwd = 1.6)
  for (nm in names(res)) {
    r <- res[[nm]]
    graphics::points(pmax(r$outdegree, 1), pmax(r$n_sig, 1), pch = 16, cex = 0.32,
                     col = cols[[nm]])
  }
  graphics::points(pmax(aut_od, 1), pmax(aut_n, 1), pch = 17, cex = 1.7, col = "black")
  ## Short label: the bound line passes through the space left of the point, so
  ## a longer label would sit on the line. The network size is in the legend of
  ## the main text and in docs/06_core_finding_CN.md.
  graphics::text(pmax(aut_od, 1), pmax(aut_n, 1),
                 "published\nTrem2", pos = 1, cex = BIB_CEX_TXT)
  ## The bound line leaves the panel at x = 3000, so the region right of it is
  ## free; the label sits there rather than on the line itself.
  bp <- usr_frac(0.97, 0.52)
  graphics::text(bp[1], bp[2], "|DR| <= outdegree + 1", cex = BIB_CEX_TXT,
                 col = "grey25", adj = 1)
  graphics::legend("bottomright", legend = bib_short(names(res)), pch = 16, cex = BIB_CEX_LEG,
                   col = cols, bty = "o", bg = "white", box.col = "white",
                   inset = 0.01, pt.cex = 0.9)

  ## Panel B
  graphics::plot(NA, xlim = c(1, 5000), ylim = c(0, 1.05), log = "x",
                 xlab = xlab_ko,
                 ylab = "Fraction of targets significant",
                 main = "B  Weight dilution", xaxt = "n")
  graphics::axis(1, at = c(1, 10, 100, 1000, 5000), labels = c("1", "10", "100", "1000", "5000"))
  for (nm in names(res)) {
    r <- res[[nm]][!is.na(res[[nm]]$frac_targets_sig), ]
    graphics::points(pmax(r$outdegree, 1), r$frac_targets_sig, pch = 16, cex = 0.32,
                     col = cols[[nm]])
  }
  graphics::points(pmax(aut_od, 1), aut_n / aut_od, pch = 17, cex = 1.7, col = "black")
  ## Outside the panel: both corners are taken (the top-right by the single
  ## published point, the bottom-left by the mass of low-outdegree points).
  graphics::mtext("low outdegree: all targets significant", side = 3, line = 0.1,
                  cex = BIB_CEX_TXT, col = "grey25", adj = 0.5)

  ## Panel C. Every bar sits at ~100%, so the value labels go into the top
  ## margin (they cannot be clipped there) and each knockout count joins the
  ## network name on the x axis, where the label already wraps to two lines.
  comp <- do.call(rbind, lapply(res, function(r) c(sum(r$all_inside), nrow(r))))
  compl <- 100 * comp[, 1] / comp[, 2]
  ## One short line per bar: a two-line label rotated by 90 degrees collided with
  ## its neighbours, and a one-line label carrying the knockout count needs more
  ## height than the axis has. The counts stay in the caption.
  net_lab <- bib_short(names(res))
  b <- graphics::barplot(compl, ylim = c(0, 112),
                         col = cols,
                         names.arg = net_lab,
                         border = NA, las = 2, cex.names = BIB_CEX_AXIS,
                         ylab = "DR genes contained (%)",
                         main = "C  Containment holds")
  graphics::abline(h = 100, lty = 2, col = "grey35")
  bib_bar_labels_rot(b, sprintf("%.1f%%", compl), line = 0.3)
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_degree_artifact"),
           width_mm = BIB_FULL_MM, height_mm = 100)
cat("\nSaved results/figures/Fig_degree_artifact.{pdf,tiff,png}\n")
