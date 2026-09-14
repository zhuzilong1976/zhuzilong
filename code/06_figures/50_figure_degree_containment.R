# Publication-quality figure: the DR gene set is determined by the KO gene's outdegree
# All labels in English to avoid font-rendering issues on Windows
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"
ndir <- "outputs/results/technical_null"
fig_dir <- "outputs/figures"
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

aut <- readRDS("work/author/author_validation.rds")
aut_od <- as.numeric(aut$od[aut$gKO])
aut_n <- aut$n_sig

draw <- function() {
  op <- graphics::par(mfrow = c(1, 3), mar = c(4.6, 4.6, 3.2, 1), oma = c(0, 0, 2.6, 0),
                      family = "sans")

  ## Panel A
  graphics::plot(NA, xlim = c(1, 20000), ylim = c(0.8, 3000), log = "xy",
                 xlab = "Outdegree of knocked-out gene",
                 ylab = "Differentially regulated genes (|DR|)",
                 main = "A  |DR| is governed by outdegree", cex.lab = 1.05, cex.main = 1.1)
  graphics::abline(0, 1, lty = 2, col = "grey35", lwd = 1.6)
  for (nm in names(res)) {
    r <- res[[nm]]
    graphics::points(pmax(r$outdegree, 1), pmax(r$n_sig, 1), pch = 16, cex = 0.32,
                     col = grDevices::adjustcolor(r$col[1], 0.4))
  }
  graphics::points(pmax(aut_od, 1), pmax(aut_n, 1), pch = 17, cex = 1.7, col = "black")
  graphics::text(pmax(aut_od, 1), pmax(aut_n, 1),
                 " published Trem2\n (7,715-gene network)", pos = 2, cex = 0.68)
  graphics::text(3.2, 2600, "|DR| <= outdegree + 1", cex = 0.72, col = "grey25", adj = 0)
  graphics::legend("bottomright", legend = names(res), pch = 16, cex = 0.6,
                   col = vapply(res, function(x) x$col[1], character(1)),
                   bty = "n", pt.cex = 0.85)

  ## Panel B
  graphics::plot(NA, xlim = c(1, 5000), ylim = c(0, 1.05), log = "x",
                 xlab = "Outdegree of knocked-out gene",
                 ylab = "Fraction of targets significant",
                 main = "B  Weight-dilution effect", cex.lab = 1.05, cex.main = 1.1, xaxt = "n")
  graphics::axis(1, at = c(1, 10, 100, 1000, 5000), labels = c("1", "10", "100", "1000", "5000"))
  for (nm in names(res)) {
    r <- res[[nm]][!is.na(res[[nm]]$frac_targets_sig), ]
    graphics::points(pmax(r$outdegree, 1), r$frac_targets_sig, pch = 16, cex = 0.32,
                     col = grDevices::adjustcolor(r$col[1], 0.4))
  }
  graphics::points(pmax(aut_od, 1), aut_n / aut_od, pch = 17, cex = 1.7, col = "black")
  graphics::text(1.6, 0.98, "low outdegree:\nconcentrated weights\n-> all targets significant",
                 cex = 0.62, adj = 0, col = "grey25")

  ## Panel C
  comp <- do.call(rbind, lapply(res, function(r) c(sum(r$all_inside), nrow(r))))
  compl <- 100 * comp[, 1] / comp[, 2]
  b <- graphics::barplot(compl, ylim = c(0, 108),
                         col = vapply(res, function(x) x$col[1], character(1)),
                         border = NA, las = 2, cex.names = 0.6,
                         ylab = "Genes with DR in {KO} U targets (%)",
                         main = "C  Containment holds almost exactly", cex.main = 1.1)
  graphics::abline(h = 100, lty = 2, col = "grey35")
  graphics::text(b, compl + 3.5, labels = sprintf("%.1f%%\n(n=%d)", compl, comp[, 2]), cex = 0.68)
  graphics::par(op)
}

png(file.path(fig_dir, "Fig_degree_artifact.png"), width = 2400, height = 900, res = 200)
draw(); dev.off()
pdf(file.path(fig_dir, "Fig_degree_artifact.pdf"), width = 12, height = 4.5)
draw(); dev.off()
cat("\nSaved outputs/figures/Fig_degree_artifact.{png,pdf}\n")
