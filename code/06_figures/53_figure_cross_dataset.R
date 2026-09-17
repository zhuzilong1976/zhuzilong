# Figure: cross-dataset replication of the containment rule
source("code/06_figures/_bib_figstyle.R")
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

our <- utils::read.csv("results/figures/containment_all_networks.csv", stringsAsFactors = FALSE)
rep <- utils::read.csv("outputs/results/replication/containment_across_datasets.csv",
                       stringsAsFactors = FALSE)

rep$od1 <- as.numeric(sub("\\+.*", "", rep$ko_outdegree))
rep$od2 <- ifelse(grepl("\\+", rep$ko_outdegree),
                  as.numeric(sub(".*\\+", "", rep$ko_outdegree)), NA)
rep$od_total <- rep$od1 + ifelse(is.na(rep$od2), 0, rep$od2)

draw <- function() {
  op <- bib_par(mfrow = c(1, 2), mar = c(5.4, 5.0, 3.0, 1.0), oma = c(0, 0, 0, 0))

  graphics::plot(NA, xlim = c(1, 12000), ylim = c(1, 1200), log = "xy",
                 xlab = "Outdegree of the knocked-out gene in the inferred network",
                 ylab = "Differentially regulated genes (|DR|)",
                 main = "A  Rule holds across 14 networks")
  graphics::abline(0, 1, lty = 2, col = "grey45", lwd = 1.6)
  graphics::points(pmax(our$outdegree, 1), pmax(our$n_sig, 1), pch = 16, cex = 0.40,
                   col = "#b0b6b8")
  graphics::points(rep$od_total, rep$n_DR, pch = 21, cex = 1.5, lwd = 1.1,
                   bg = "#e74c3c", col = "black")
  ## Points with the same gene symbol (Mecp2, Cftr) are labelled on opposite
  ## sides so the labels do not collide.
  lab_pos <- rep(4, nrow(rep))
  for (g in unique(rep$ko_genes[duplicated(rep$ko_genes) |
                                duplicated(rep$ko_genes, fromLast = TRUE)])) {
    idx <- which(rep$ko_genes == g)
    lab_pos[idx] <- c(4, 2)[seq_along(idx) %% 2 + 1]
  }
  ## Points with a nearly identical outdegree (Dmd, Malat1, Nkx2-1) would put
  ## their labels on top of each other; cycle below / above / left instead.
  ord <- order(rep$od_total)
  cl <- cumsum(c(TRUE, diff(log10(rep$od_total[ord])) > 0.15))
  for (c_id in unique(cl)) {
    idx <- ord[cl == c_id]
    if (length(idx) > 1) lab_pos[idx] <- c(4, 1, 3, 2)[seq_along(idx) %% 4 + 1]
  }
  graphics::text(rep$od_total, rep$n_DR, labels = rep$ko_genes, pos = lab_pos,
                 cex = BIB_CEX_TXT, offset = 0.4)
  graphics::text(1.4, 700, "|DR| <= outdegree + 1", adj = 0, cex = BIB_CEX_TXT,
                 col = "grey25")
  graphics::legend("bottomright",
                   legend = c("this study\n(3,195 knockouts)",
                              "published\n(10 datasets)"),
                   pch = c(16, 21), pt.bg = c(NA, "#e74c3c"), col = c("#7f8c8d", "black"),
                   pt.cex = c(1, 1.3), bty = "o", bg = "white", box.col = "white",
                   cex = BIB_CEX_LEG)

  our_ok <- our$all_inside %in% c(TRUE, "TRUE", "true")
  net_tab <- tapply(our_ok, our$network, mean)
  ## Gene symbols only: the datasets behind each bar are listed in the
  ## manuscript figure legend and in Table S3.
  labs <- c(bib_short(names(net_tab)), rep$ko_genes)
  pct <- c(100 * as.numeric(net_tab), rep$pct_inside)
  cat("names(net_tab):", length(net_tab), "| rep rows:", nrow(rep),
      "| length(pct):", length(pct), "| length(labs):", length(labs), "\n")
  grp <- c(rep("this study", nrow(our)), rep("published", nrow(rep)))
  cols <- ifelse(grp == "this study", "#7f8c8d", "#e74c3c")
  b <- graphics::barplot(pct, names.arg = labs, col = cols, border = NA, ylim = c(0, 128),
                         las = 2, cex.names = BIB_CEX_AXIS,
                         ylab = "DR genes inside {KO} U direct targets (%)",
                         main = "B  Containment: 100% in every dataset")
  graphics::abline(h = 100, lty = 2, col = "grey40")
  graphics::text(b, pct + 3, labels = sprintf("%.1f", pct), cex = BIB_CEX_TXT, srt = 90)
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_replication"),
           width_mm = BIB_FULL_MM, height_mm = 100)
cat("Saved results/figures/Fig_replication.{pdf,tiff,png}\n")
