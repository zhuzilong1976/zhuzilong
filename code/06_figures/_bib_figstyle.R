# =============================================================================
# Shared output settings for every submission figure.
#
# Target: OUP / Briefings in Bioinformatics figure requirements
#   * full-page-width artwork, 170 mm (single column would be 85 mm)
#   * sans-serif font, minimum 8 pt at final printed size
#   * fonts embedded in the PDF
#   * RGB raster backups at >= 300 dpi (600 dpi for combined line/halftone)
#   * one figure per file, no transparency, white background
#
# Usage inside a figure script:
#   source("code/06_figures/_bib_figstyle.R")
#   draw <- function() { ... }          # unchanged drawing code
#   bib_render(draw, file.path(fig_dir, "Fig_degree_artifact"),
#              width_mm = BIB_FULL_MM, height_mm = 72)
#
# Every device is opened with the same point size and font, so the PDF, the
# TIFF and the PNG show identical typography.
# =============================================================================

BIB_FULL_MM <- 170        # full page width used by the journal
BIB_SINGLE_MM <- 85       # single column width, available if a figure is split

# Base point size. Every cex in the drawing code must be >= 0.9, so the
# smallest text in a figure is 0.9 * 9 = 8.1 pt, i.e. above the 8 pt floor.
BIB_PS <- 9

BIB_FONT <- "Arial"       # embedded in the PDF, matches Helvetica metrics

## Standard cex values -------------------------------------------------------
## Values below BIB_CEX_MIN are not allowed: they would print under 8 pt.
BIB_CEX_MIN  <- 0.9
BIB_CEX_MAIN <- 1.15      # 10.35 pt panel titles
BIB_CEX_LAB  <- 1.00      # 9.0 pt axis titles
BIB_CEX_AXIS <- 0.90      # 8.1 pt tick labels
BIB_CEX_TXT  <- 0.90      # 8.1 pt in-panel annotations
BIB_CEX_LEG  <- 0.90      # 8.1 pt legends

## Solid tints that replace grDevices::adjustcolor(..., alpha) so that no
## transparency reaches the typesetter.
BIB_TINT <- list(
  blue  = "#8fbcdc", lightblue = "#bcd8ec",
  red   = "#e0a49c", orange    = "#f0c8a0",
  grey  = "#b8b8b8", darkgrey  = "#8c8c8c",
  purple = "#c3a8cf", teal     = "#a7d3cc"
)

bib_par <- function(...) {
  ## Setting mfrow/mfcol makes R shrink par("cex") below 1 (0.66 for a
  ## three-panel row), which would print every label under 6 pt. The layout is
  ## therefore applied first and the style second, with cex pinned to 1.
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(...)
  graphics::par(family = BIB_FONT, bg = "white", fg = "black",
                cex = 1,
                col.axis = "black", col.lab = "black", col.main = "black",
                cex.axis = BIB_CEX_AXIS, cex.lab = BIB_CEX_LAB,
                cex.main = BIB_CEX_MAIN)
  invisible(op)
}

## Display names for the four networks of this study. Full names stay in
## results/figures/containment_all_networks.csv and in the manuscript legend;
## the panels use the short forms so that rotated labels stay inside the panel.
BIB_NET_SHORT <- c("MS microglia, control net" = "Control net",
                   "MS microglia, lesion net"   = "MS lesion net",
                   "Control split-half A"       = "Split-half A",
                   "Control split-half B"       = "Split-half B")

bib_short <- function(x) {
  out <- unname(BIB_NET_SHORT[x])
  ifelse(is.na(out), x, out)
}

## Open the three output devices with identical geometry and typography.
.bib_open <- function(stem, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  list(
    pdf = function() grDevices::cairo_pdf(
      paste0(stem, ".pdf"), width = w, height = h, pointsize = BIB_PS,
      family = BIB_FONT, bg = "white", onefile = FALSE),
    tiff = function() grDevices::tiff(
      paste0(stem, ".tiff"), width = w, height = h, units = "in",
      pointsize = BIB_PS, res = 600, compression = "lzw",
      family = BIB_FONT, bg = "white", type = "cairo"),
    png = function() grDevices::png(
      paste0(stem, ".png"), width = w, height = h, units = "in",
      pointsize = BIB_PS, res = 300,
      family = BIB_FONT, bg = "white", type = "cairo")
  )
}

## Draw the same function to PDF (submission copy), TIFF (600 dpi raster
## backup) and PNG (300 dpi preview). Returns the written paths.
bib_render <- function(draw, stem, width_mm = BIB_FULL_MM, height_mm) {
  stopifnot(is.function(draw), length(height_mm) == 1)
  dev <- .bib_open(stem, width_mm, height_mm)
  out <- character(0)
  for (nm in names(dev)) {
    dev[[nm]]()
    op <- graphics::par(no.readonly = TRUE)   # restore par for every device
    draw()
    graphics::par(op)
    grDevices::dev.off()
    out <- c(out, paste0(stem, ".", nm))
  }
  cat(sprintf("  %s  (%.0f x %.0f mm, base %g pt, min text %g pt)\n",
              basename(stem), width_mm, height_mm, BIB_PS,
              BIB_PS * BIB_CEX_MIN))
  invisible(out)
}
