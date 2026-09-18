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
  ## Right margin 1.6, not 1.0: panel A's y-axis title is a physical width, and
  ## at 1.0 the last character landed outside the device and was cut off.
  ## mar bottom 6.6: panel B's two-line rotated dataset labels ran past the
  ## bottom edge at 5.4 (ink on the last raster row). mar right 2.0 also keeps
  ## panel B's long x-axis title inside the device.
  op <- bib_par(mfrow = c(1, 2), mar = c(8.6, 5.0, 3.0, 2.0), oma = c(0, 0, 0, 0))

  ## The y axis runs to 1,200 although the highest point is 377: the extra
  ## decade is headroom for labels. Nkx2-1, Dmd and Malat1 lie within 0.06 of a
  ## decade of each other on the y axis and within 0.21 on x, so their three
  ## labels cannot be fitted around the points at all; they are stacked in the
  ## empty band above the cluster instead.
  graphics::plot(NA, xlim = c(1, 20000), ylim = c(1, 1200), log = "xy",
                 xlab = "Outdegree in the inferred network",
                 ylab = "absolute DR genes",
                 main = "A  Rule holds across 14 networks")
  graphics::abline(0, 1, lty = 2, col = "grey45", lwd = 1.6)
  graphics::points(pmax(our$outdegree, 1), pmax(our$n_sig, 1), pch = 16, cex = 0.40,
                   col = "#b0b6b8")
  graphics::points(rep$od_total, rep$n_DR, pch = 21, cex = 1.5, lwd = 1.1,
                   bg = "#e74c3c", col = "black")
  ## Label centres, in the order of rep$ko_genes: Hnf4a+Hnf4g, Nkx2-1, Dmd,
  ## Malat1, Ahr, Mecp2, Mecp2, Cftr, Akap7, Cftr.
  ##
  ## At 8.1 pt a label is about 1.0 log unit wide and 0.13 log units tall, so two
  ## labels collide unless they differ by more than that on at least one axis.
  ## These ten positions were checked one against another under that rule. They
  ## replace a set in which Hnf4a+Hnf4g sat on top of the bound-line note and the
  ## two Mecp2 labels plus Nkx2-1 ran together at the top of the cluster.
  lab_x <- c( 126, 8912, 4000, 3548,  224, 1100,  330,   56,  4.2, 6.3)
  lab_y <- c(  96,  263,  700,  105,   40,  230,  140,   16,  2.6,  24)
  graphics::text(lab_x, lab_y, labels = rep$ko_genes, cex = BIB_CEX_TXT)
  ## One legend for the whole panel, in the empty upper-left triangle, one short
  ## line per entry. The bound line has its own entry, so no note has to sit in
  ## the plotting region, where the old one ran into the Hnf4a+Hnf4g label. The
  ## dataset counts stay in the figure legend: spelling them out here made the box
  ## 77% of the panel wide and it covered the Mecp2 and Dmd labels.
  graphics::legend("topleft",
                   legend = c("this study",
                              "published",
                              paste0("|DR| ", "\u2264", " outdegree + 1")),
                   pch = c(16, 21, NA), lty = c(0, 0, 2),
                   pt.bg = c(NA, "#e74c3c", NA), col = c("#7f8c8d", "black", "grey45"),
                   pt.cex = c(1, 1.3, NA), lwd = c(1, 1, 1.6),
                   bty = "o", bg = "white", box.col = "white",
                   cex = BIB_CEX_LEG)

  our_ok <- our$all_inside %in% c(TRUE, "TRUE", "true")
  net_tab <- tapply(our_ok, our$network, mean)
  n_net <- tapply(our_ok, our$network, length)
  ## tapply sorts the names alphabetically, which would put the two split-halves
  ## before the two networks they are split from. Keep the order used in the
  ## manuscript: Control, MS, split-half A, split-half B.
  net_order <- c("MS microglia, control net", "MS microglia, lesion net",
                 "Control split-half A", "Control split-half B")
  net_tab <- net_tab[net_order]; n_net <- n_net[net_order]
  ## Network size joins the axis label, which already wraps to two lines; the
  ## datasets behind each bar stay in the manuscript legend and Table S3.
  ## One short line per label: fourteen two-line rotated labels in a 74 mm panel
  ## overlap each other. Network sizes stay in the figure legend and Table S3.
  labs <- c(bib_short(names(net_tab)), rep$ko_genes)
  pct <- c(100 * as.numeric(net_tab), rep$pct_inside)
  ## grp has to have one entry per bar. Spelled with nrow(our) it had 3,195
  ## entries, so the first 14 were all "this study" and every bar came out grey,
  ## including the 10 published datasets that are meant to be red.
  grp <- c(rep("this study", length(net_tab)), rep("published", nrow(rep)))
  cols <- ifelse(grp == "this study", "#7f8c8d", "#e74c3c")
  b_rep <- graphics::barplot(pct, names.arg = rep("", length(labs)), col = cols,
                             border = NA, ylim = c(0, 112), las = 1,
                             ylab = "DR genes inside {KO} U targets (%)",
                             main = "B  Containment: 100% in every dataset")
  graphics::abline(h = 100, lty = 2, col = "grey40")
  ## Dataset labels read from the top down, as in Figures 1C and S1, instead of
  ## R's bottom-up default.
  usr <- graphics::par("usr")
  graphics::text(b_rep, usr[3] - 0.015 * (usr[4] - usr[3]), labels = labs,
                 srt = -90, adj = c(0, 0.5), cex = BIB_CEX_AXIS, xpd = NA)
  ## The bars do not repeat 14 identical values above the line: the exact rate
  ## is in the panel title, the figure legend and Table S3.
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_replication"),
           width_mm = BIB_FULL_MM, height_mm = 110)
cat("Saved results/figures/Fig_replication.{pdf,tiff,png}\n")
