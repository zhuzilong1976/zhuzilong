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
## ---------------------------------------------------------------------------
## Headroom for labels drawn above a bar.
##
## A bar chart with a fixed ylim and labels above the bar clips them as soon as
## the bars approach the top of the axis, which is exactly what happens in the
## containment figures where every bar sits at ~100%.
##
## A label rotated by srt = 90 and anchored with adj = c(0, 0.5) starts at the
## given y and runs upward for its own width. The conversion needs care: the
## value of strwidth(units = "user") is expressed in x-axis user units, so the
## height of one line of text in *y* units has to be measured first and used as
## the in/inch factor. par("pin") gives the plot height in inches.
usr_rot_label_h <- function(labels, cex = BIB_CEX_TXT) {
  usr <- graphics::par("usr")
  din <- graphics::par("din")           # device size in inches
  if (din[2] <= 0) return(0)            # no device geometry
  old <- graphics::par(cex = cex)
  on.exit(graphics::par(old), add = TRUE)
  ## The label turns through 90 degrees, so the space it needs along y equals
  ## its own printed width. strwidth("X") gives the width of one character in
  ## the same units, so the ratio is the label width in character advances;
  ## strheight("X") converts a character advance into y user units.
  max(graphics::strwidth(labels)) / graphics::strwidth("X") *
    graphics::strheight("X")
}

## Value labels for a bar chart whose bars all sit near the top of the axis.
##
## Drawing them as text() above the bar needs the axis to reserve headroom that
## has to be measured before the panel exists, which is unreliable. mtext()
## writes into the figure margin instead, so the labels can never be clipped
## and the axis range stays exactly as specified. `at` is in the same user
## coordinates as the barplot() return value.
bib_bar_labels <- function(x, labels, line = 0.15, cex = BIB_CEX_TXT,
                           col = "black", font = 1, srt = 0) {
  for (i in seq_along(x)) {
    graphics::mtext(labels[i], side = 3, at = x[i], line = line,
                    cex = cex, col = col, font = font, xpd = NA, srt = srt)
  }
}
## Alias: the same margin trick, but the labels run vertically, which is what a
## bar chart needs when the bars are too narrow to fit the values side by side.
bib_bar_labels_rot <- function(x, labels, line = 0.15, cex = BIB_CEX_TXT,
                               col = "black", font = 1) {
  bib_bar_labels(x, labels, line = line, cex = cex, col = col, font = font,
                 srt = 90)
}

## Height of `n_lines` lines of text in y-axis user units (for placing a value
## label inside a bar, or a multi-line annotation at a computed height).
usr_text_h <- function(n_lines = 1, cex = BIB_CEX_TXT) {
  old <- graphics::par(cex = cex)
  on.exit(graphics::par(old), add = TRUE)
  n_lines * graphics::strheight("X", units = "user")
}

## Screen-space (normalised device coordinate) extent of a label drawn at (x, y)
## with the given justification. Used to test two labels for overlap without
## needing the bounding-box API that base graphics does not expose.
## Place a label next to each point, choosing the first compass direction that
## does not overlap an already-placed label or a data point. Cardinal
## directions are tried before diagonals, so labels stay axis-aligned wherever
## the figure allows it.
##
## This exists because fixed direction rules (below, above, left, repeat) break
## when several points sit within a character of each other on both axes, which
## is exactly the case for the Dmd / Malat1 / Nkx2-1 cluster in Figure 2A.
##
## Geometry is done in normalised device coordinates. The conversion avoids the
## two traps in base graphics: strwidth(units = "user") is an *x* length (so it
## is wrong on a log or differently scaled y axis), and par("din") is in inches
## while strwidth() is in user units. Both widths are therefore expressed as
## multiples of one character advance, and a single measured y advance converts
## to inches, which gives the x scale as well.
bib_point_labels <- function(x, y, labels, cex = BIB_CEX_TXT,
                             offset_in = 0.03, avoid_points = TRUE,
                             point_cex = 1.5) {
  stopifnot(length(x) == length(y), length(x) == length(labels))
  lim <- graphics::par("usr")
  din <- graphics::par("din")
  old <- graphics::par(cex = cex)
  on.exit(graphics::par(old), add = TRUE)

  ## On a log axis par("usr") and the placement coordinates are log10 values
  ## while x and y are raw data, so everything is done in the par("usr") space
  ## and converted back once at the end.
  log_x <- graphics::par("xlog")
  log_y <- graphics::par("ylog")
  ux <- if (log_x) log10(x) else x
  uy <- if (log_y) log10(y) else y

  ## Label size in inches. Reading it off the device is the only reliable way:
  ## strwidth(units = "user") is an x-axis length and is wrong on a log axis,
  ## and par("din") is in inches while strwidth() is not.
  char_h_in <- graphics::par("cin")[2]            # one character, inches
  text_w_in <- vapply(labels, function(s) {
    graphics::strwidth(s) / graphics::strwidth("X") * char_h_in
  }, numeric(1))
  plot_w_in <- din[1] - sum(graphics::par("mai")[c(2, 4)])
  plot_h_in <- din[2] - sum(graphics::par("mai")[c(1, 3)])

  per_in_x <- (lim[2] - lim[1]) / plot_w_in       # usr units per inch
  per_in_y <- (lim[4] - lim[3]) / plot_h_in
  w_u <- text_w_in * per_in_x                     # label width, usr units
  h_u <- char_h_in * per_in_y
  pt_u <- graphics::par("cin")[1] * point_cex / 2 * per_in_x
  pt_v <- graphics::par("cin")[1] * point_cex / 2 * per_in_y
  off <- offset_in * per_in_x

  dirs <- list(c(0, 1), c(0, -1), c(1, 0), c(-1, 0),
               c(1, 1), c(-1, 1), c(1, -1), c(-1, -1))
  boxes <- vector("list", length(x))
  out_x <- ux; out_y <- uy

  overlaps <- function(b) {
    for (j in seq_along(boxes)) {
      q <- boxes[[j]]
      if (!is.null(q) && b[1] < q[2] && b[2] > q[1] &&
          b[3] < q[4] && b[4] > q[3]) return(TRUE)
    }
    if (avoid_points) {
      for (k in seq_along(ux)) {
        if (b[1] < ux[k] + pt_u && b[2] > ux[k] - pt_u &&
            b[3] < uy[k] + pt_v && b[4] > uy[k] - pt_v) return(TRUE)
      }
    }
    FALSE
  }

  for (i in seq_along(ux)) {
    w <- w_u[i] / 2
    h <- h_u / 2
    placed <- NULL
    for (d in dirs) {
      cx <- ux[i] + d[1] * (w + off)
      cy <- uy[i] + d[2] * (h + off)
      b <- c(cx - w, cx + w, cy - h, cy + h)
      if (b[1] < lim[1] || b[2] > lim[2] ||
          b[3] < lim[3] || b[4] > lim[4]) next
      if (overlaps(b)) next
      placed <- list(cx = cx, cy = cy, box = b)
      break
    }
    if (is.null(placed)) {
      cx <- ux[i]
      cy <- uy[i] + h + off
      placed <- list(cx = cx, cy = cy, box = c(cx - w, cx + w, cy - h, cy + h))
    }
    boxes[[i]] <- placed$box
    out_x[i] <- placed$cx
    out_y[i] <- placed$cy
  }
  list(x = if (log_x) 10^out_x else out_x,
       y = if (log_y) 10^out_y else out_y)
}

## A point given as a fraction of the current plot region, so an annotation
## stays inside the panel when the axis range changes.
usr_frac <- function(x_frac, y_frac) {
  usr <- graphics::par("usr")
  c(x_frac * (usr[2] - usr[1]) + usr[1], y_frac * (usr[4] - usr[3]) + usr[3])
}

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
