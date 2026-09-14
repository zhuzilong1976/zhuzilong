# Figure: cross-dataset replication of the containment rule
fig_dir <- "outputs/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

our <- utils::read.csv("outputs/figures/containment_all_networks.csv", stringsAsFactors = FALSE)
rep <- utils::read.csv("outputs/results/replication/containment_across_datasets.csv",
                       stringsAsFactors = FALSE)

rep$od1 <- as.numeric(sub("\\+.*", "", rep$ko_outdegree))
rep$od2 <- ifelse(grepl("\\+", rep$ko_outdegree),
                  as.numeric(sub(".*\\+", "", rep$ko_outdegree)), NA)
rep$od_total <- rep$od1 + ifelse(is.na(rep$od2), 0, rep$od2)

draw <- function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.2, 1), oma = c(0, 0, 2.2, 0),
                      family = "sans")

  graphics::plot(NA, xlim = c(1, 12000), ylim = c(1, 1200), log = "xy",
                 xlab = "Outdegree of the knocked-out gene in the inferred network",
                 ylab = "Differentially regulated genes (|DR|)",
                 main = "A  Rule holds across 14 networks", cex.main = 1.05, cex.lab = 0.95)
  graphics::abline(0, 1, lty = 2, col = "grey45", lwd = 1.6)
  graphics::points(pmax(our$outdegree, 1), pmax(our$n_sig, 1), pch = 16, cex = 0.28,
                   col = grDevices::adjustcolor("#7f8c8d", 0.30))
  graphics::points(rep$od_total, rep$n_DR, pch = 21, cex = 1.5, lwd = 1.1,
                   bg = "#e74c3c", col = "black")
  graphics::text(rep$od_total, rep$n_DR, labels = rep$ko_genes, pos = 4, cex = 0.62, offset = 0.4)
  graphics::text(1.4, 700, "|DR| <= outdegree + 1", adj = 0, cex = 0.8, col = "grey25")
  graphics::legend("bottomright",
                   legend = c("this study: 3,195 KOs in 4 networks (MS microglia)",
                              "published by Cai lab: 10 independent datasets"),
                   pch = c(16, 21), pt.bg = c(NA, "#e74c3c"), col = c("#7f8c8d", "black"),
                   pt.cex = c(1, 1.3), bty = "n", cex = 0.68)

  our_ok <- our$all_inside %in% c(TRUE, "TRUE", "true")
  net_tab <- tapply(our_ok, our$network, mean)
  labs <- c(names(net_tab), paste0(rep$ko_genes, " | ", sub(" / .*", "", rep$dataset)))
  pct <- c(100 * as.numeric(net_tab), rep$pct_inside)
  cat("names(net_tab):", length(net_tab), "| rep rows:", nrow(rep),
      "| length(pct):", length(pct), "| length(labs):", length(labs), "\n")
  grp <- c(rep("this study", nrow(our)), rep("published", nrow(rep)))
  cols <- ifelse(grp == "this study", "#7f8c8d", "#e74c3c")
  b <- graphics::barplot(pct, names.arg = labs, col = cols, border = NA, ylim = c(0, 112),
                         las = 2, cex.names = 0.55,
                         ylab = "DR genes inside {KO} U direct targets (%)",
                         main = "B  Containment: 100% in every dataset", cex.main = 1.05)
  graphics::abline(h = 100, lty = 2, col = "grey40")
  graphics::text(b, pct + 3, labels = sprintf("%.1f", pct), cex = 0.6, srt = 90)
  graphics::par(op)
}

png(file.path(fig_dir, "Fig_replication.png"), width = 2400, height = 1000, res = 190)
draw(); dev.off()
pdf(file.path(fig_dir, "Fig_replication.pdf"), width = 12.5, height = 5.2)
draw(); dev.off()
cat("Saved outputs/figures/Fig_replication.{png,pdf}\n")
