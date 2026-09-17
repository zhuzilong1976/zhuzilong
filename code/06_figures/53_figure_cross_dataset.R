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
  op <- bib_par(mfrow = c(1, 2), mar = c(5.4, 5.0, 3.0, 1.6), oma = c(0, 0, 0, 0))

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
  ## Label positions, in the order of rep$ko_genes: Hnf4a+Hnf4g, Nkx2-1, Dmd,
  ## Malat1, Ahr, Mecp2, Mecp2, Cftr, Akap7, Cftr. Columns are the label's
  ## x and y in data units; both are given explicitly because no rule fitted
  ## them, as the note on the y range above explains.
  lab_x <- c(  60,  3600, 2000, 2900,  60,   700,  100,  38,    2.2,  6)
  lab_y <- c(  27,   700,  500,  340,  30,   640,  470,  17,    3.2, 11)
  graphics::text(lab_x, lab_y, labels = rep$ko_genes, cex = BIB_CEX_TXT)
  ## Above the bound line at the far left: the line is still below y = 90 where
  ## the panel starts, and no point is plotted in that corner.
  graphics::text(1.15, 120, "|DR| <= outdegree + 1", adj = 0, cex = BIB_CEX_TXT,
                 col = "grey25")
  graphics::legend("bottomright",
                   legend = c("this study\n(3,195 knockouts)",
                              "published\n(10 datasets)"),
                   pch = c(16, 21), pt.bg = c(NA, "#e74c3c"), col = c("#7f8c8d", "black"),
                   pt.cex = c(1, 1.3), bty = "o", bg = "white", box.col = "white",
                   cex = BIB_CEX_LEG)

  our_ok <- our$all_inside %in% c(TRUE, "TRUE", "true")
  net_tab <- tapply(our_ok, our$network, mean)
  n_net <- tapply(our_ok, our$network, length)
  ## Network size joins the axis label, which already wraps to two lines; the
  ## datasets behind each bar stay in the manuscript legend and Table S3.
  labs <- c(sprintf("%s\n(n=%d)", bib_short(names(net_tab)), n_net),
            sprintf("%s\n(n=%s)", rep$ko_genes, gsub("[+]", "+\n", rep$network_genes)))
  pct <- c(100 * as.numeric(net_tab), rep$pct_inside)
  grp <- c(rep("this study", nrow(our)), rep("published", nrow(rep)))
  cols <- ifelse(grp == "this study", "#7f8c8d", "#e74c3c")
  graphics::barplot(pct, names.arg = labs, col = cols, border = NA,
                    ylim = c(0, 112), las = 2, cex.names = BIB_CEX_AXIS,
                    ylab = "DR genes inside {KO} U targets (%)",
                    main = "B  Containment: 100% in every dataset")
  graphics::abline(h = 100, lty = 2, col = "grey40")
  ## The bars do not repeat 14 identical values above the line: the exact rate
  ## is in the panel title, the figure legend and Table S3.
  graphics::par(op)
}

bib_render(draw, file.path(fig_dir, "Fig_replication"),
           width_mm = BIB_FULL_MM, height_mm = 100)
cat("Saved results/figures/Fig_replication.{pdf,tiff,png}\n")
