# ============================================================
# split_dual_tissue.R
#
# For Visium samples with TWO tissue sections placed on one capture
# area (e.g. GSM9622958 / "0223", pre+post treatment on one slide):
# clusters spots by physical position to separate the two tissues,
# then exports EACH piece as its own independent Space-Ranger-style
# folder - so your downstream Python script (cell2location etc.) can
# process them as two completely separate samples with zero changes.
#
# Uses dbscan (density-based clustering) rather than a manually built
# neighbor graph: it naturally finds dense groups of spots separated
# by empty background, and cleanly flags stray/edge spots as "noise"
# instead of forcing them into either tissue piece.
# ============================================================

# Locate this script's own folder so the source() call below works regardless
# of your current working directory (fixes "cannot open the connection" if
# you're sourcing this from outside the visium_scripts/ folder).
.split_dual_tissue_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NULL
)
if (is.null(.split_dual_tissue_dir) || !nzchar(.split_dual_tissue_dir)) {
  # Fallback if sourced in a way that doesn't expose the file path (e.g. pasted
  # directly into the console) - assumes working directory IS visium_scripts/
  .split_dual_tissue_dir <- "."
}
source(file.path(.split_dual_tissue_dir, "export_rda_to_spaceranger.R"))  # brings in write_seurat_to_spaceranger() and get_full_res_coords()
library(Seurat)

if (!requireNamespace("dbscan", quietly = TRUE)) {
  stop("Package 'dbscan' is required for spatial clustering.\n",
       "Install it with: install.packages('dbscan')")
}

#' Cluster spots by physical position to identify separate tissue pieces.
#' Does NOT export anything - just labels + lets you visually check first.
#'
#' @param obj Seurat spatial object containing both tissue pieces
#' @param image_name which image slot to use (defaults to the first one)
#' @param eps dbscan neighborhood radius, in the same pixel units as
#'   imagerow/imagecol. If NULL, auto-estimated from the data (recommended
#'   to start with the auto value, then adjust after checking the plot).
#' @param min_pts dbscan minPts - minimum spots to seed a cluster (dbscan default-ish)
#' @param min_component_size clusters smaller than this are treated as noise/debris,
#'   not a real tissue piece (guards against small artifacts getting exported as
#'   a "sample")
#' @return a list: $obj (original object with a new "tissue_piece" metadata
#'   column, cluster 0 = noise/unclustered), $cluster_sizes (table for a quick look)
label_tissue_pieces <- function(obj, image_name = NULL, eps = NULL,
                                min_pts = 5, min_component_size = 50) {
  
  if (is.null(image_name)) image_name <- Images(obj)[1]
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  xy <- as.matrix(coords[, c("imagerow", "imagecol")])
  
  if (is.null(eps)) {
    # Estimate typical spot-to-spot spacing via 6th-nearest-neighbor distance
    # (Visium spots sit on a hex grid, so ~6 immediate neighbors is expected),
    # then set eps a bit above that so adjacent spots link up but the gap
    # between the two tissues does not.
    knn_dist <- dbscan::kNNdist(xy, k = 6)
    typical_spacing <- median(knn_dist)
    eps <- typical_spacing * 2.5
    cat("Auto-estimated eps =", round(eps, 1),
        "pixels (median 6-NN spacing =", round(typical_spacing, 1), ")\n")
    cat("If this splits into the wrong number of pieces, re-run with an explicit\n",
        "'eps' argument - smaller eps = stricter separation, larger eps = more lenient.\n", sep = "")
  }
  
  clust <- dbscan::dbscan(xy, eps = eps, minPts = min_pts)
  
  # map back onto the FULL object (including any cells outside this image's coords,
  # which shouldn't normally happen but we guard for it) - default to 0 (noise)
  cluster_col <- setNames(rep(0L, ncol(obj)), colnames(obj))
  cluster_col[coords$barcode] <- clust$cluster
  obj$tissue_piece <- factor(cluster_col)
  
  sizes <- table(obj$tissue_piece)
  cat("\nCluster sizes (0 = noise/unclustered spots):\n")
  print(sizes)
  
  small_real_clusters <- names(sizes)[sizes < min_component_size & names(sizes) != "0"]
  if (length(small_real_clusters) > 0) {
    cat("\nNote: cluster(s)", paste(small_real_clusters, collapse = ", "),
        "are below min_component_size =", min_component_size,
        "- likely debris/edge artifacts rather than a real tissue piece.\n",
        "These will be excluded if you proceed to split_and_export_dual_tissue()\n")
  }
  
  list(obj = obj, cluster_sizes = sizes)
}

#' Visual sanity check BEFORE exporting - always run this first.
#' Plots spots colored by cluster over the tissue image.
plot_tissue_pieces <- function(obj, image_name = NULL) {
  if (is.null(image_name)) image_name <- Images(obj)[1]
  if (!"tissue_piece" %in% colnames(obj@meta.data)) {
    stop("Run label_tissue_pieces() first - no 'tissue_piece' column found.")
  }
  SpatialDimPlot(obj, group.by = "tissue_piece", images = image_name, pt.size.factor = 1.6) +
    ggplot2::labs(title = "Tissue piece assignment - check this looks right before exporting!")
}

#' Draw the actual split as an outline directly on the raw TIFF image -
#' more useful than colored dots for irregularly-shaped tissue, since a
#' clean boundary line makes it obvious whether the split makes sense.
#' Uses a convex hull around each cluster's spots.
#'
#' @param obj Seurat object already run through label_tissue_pieces()
#' @param tiff_path path to the ORIGINAL full-resolution TIFF (not Seurat's
#'   internal, possibly-downsampled image slot) - e.g. the raw image file
#'   you fed into Space Ranger. Coordinates in obj's coordinates slot are
#'   already in full-res pixel space, matching this file directly.
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @param max_dim if the TIFF's largest dimension exceeds this, downsample
#'   just for DISPLAY (plotting a multi-gigapixel TIFF directly is slow/
#'   memory-heavy) - point coordinates are scaled to match automatically.
#'   Set to Inf to disable and plot at full resolution.
#' @param point_size size of the individual spot markers (set to 0 to show
#'   only the hull outlines, no dots)
plot_split_boundary_on_tiff <- function(obj, tiff_path, image_name = NULL,
                                        max_dim = 2000, point_size = 0.4) {
  
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("Package 'tiff' is required to read the raw TIFF.\n",
         "Install it with: install.packages('tiff')")
  }
  if (!"tissue_piece" %in% colnames(obj@meta.data)) {
    stop("Run label_tissue_pieces() first - no 'tissue_piece' column found.")
  }
  if (!file.exists(tiff_path)) {
    stop("TIFF not found: ", tiff_path)
  }
  
  if (is.null(image_name)) image_name <- Images(obj)[1]
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  coords$tissue_piece <- as.character(obj$tissue_piece[coords$barcode])
  
  cat("Reading TIFF (this can take a moment for large full-res slide images)...\n")
  img <- tiff::readTIFF(tiff_path)
  img_h <- dim(img)[1]
  img_w <- dim(img)[2]
  cat("  image dimensions:", img_w, "x", img_h, "px\n")
  
  # Downsample just for plotting speed/memory - full-res whole-slide TIFFs
  # can be 20000+ px per side, which is far more than a screen needs to show.
  coords$plot_row <- coords$imagerow
  coords$plot_col <- coords$imagecol
  if (is.finite(max_dim) && max(img_h, img_w) > max_dim) {
    scale_factor <- max_dim / max(img_h, img_w)
    new_h <- max(1, round(img_h * scale_factor))
    new_w <- max(1, round(img_w * scale_factor))
    row_idx <- round(seq(1, img_h, length.out = new_h))
    col_idx <- round(seq(1, img_w, length.out = new_w))
    img <- img[row_idx, col_idx, , drop = FALSE]
    coords$plot_row <- coords$imagerow * scale_factor
    coords$plot_col <- coords$imagecol * scale_factor
    img_h <- new_h
    img_w <- new_w
    cat("  downsampled to", img_w, "x", img_h, "px for display (points scaled to match)\n")
  }
  
  clusters <- setdiff(unique(coords$tissue_piece), "0")
  if (length(clusters) == 0) {
    stop("No non-noise clusters found in 'tissue_piece' - nothing to draw a boundary around.")
  }
  palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
  if (length(clusters) > length(palette)) {
    palette <- grDevices::rainbow(length(clusters))
  }
  
  op <- graphics::par(mar = c(0, 0, 2, 0), mfrow = c(1, 1))
  on.exit(graphics::par(op))
  plot(0, 0, type = "n", xlim = c(0, img_w), ylim = c(img_h, 0), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = "Tissue split boundary (convex hull per cluster)")
  graphics::rasterImage(img, 0, img_h, img_w, 0)
  
  legend_labels <- character(0)
  legend_cols <- character(0)
  legend_pch <- integer(0)
  
  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    sub <- coords[coords$tissue_piece == cl, ]
    if (point_size > 0) {
      graphics::points(sub$plot_col, sub$plot_row, col = palette[i], pch = 16, cex = point_size)
    }
    if (nrow(sub) >= 3) {
      hull_idx <- grDevices::chull(sub$plot_col, sub$plot_row)
      graphics::polygon(sub$plot_col[hull_idx], sub$plot_row[hull_idx],
                        border = palette[i], lwd = 3, col = NA)
    }
    legend_labels <- c(legend_labels, paste("piece", cl, paste0("(n=", nrow(sub), ")")))
    legend_cols <- c(legend_cols, palette[i])
    legend_pch <- c(legend_pch, 16)
  }
  
  noise <- coords[coords$tissue_piece == "0", ]
  if (nrow(noise) > 0) {
    graphics::points(noise$plot_col, noise$plot_row, col = "grey40", pch = 4, cex = point_size * 1.5)
    legend_labels <- c(legend_labels, paste0("noise/excluded (n=", nrow(noise), ")"))
    legend_cols <- c(legend_cols, "grey40")
    legend_pch <- c(legend_pch, 4)
  }
  
  graphics::legend("topright", legend = legend_labels, col = legend_cols, pch = legend_pch,
                   bty = "o", bg = "white", cex = 0.8)
  
  invisible(NULL)
}

#' MANUAL split: click two points directly on the tissue image to define a
#' straight dividing line, rather than relying on automated clustering.
#' Use this when label_tissue_pieces()/dbscan can't cleanly separate the two
#' tissues (e.g. they're close together, touching, or bridged by stray spots).
#'
#' Works by testing which side of the line (point1 -> point2) each spot falls
#' on - a simple, robust approach as long as the two tissues can be separated
#' by ONE straight line (works for diagonal lines too, not just horizontal/
#' vertical cuts). See the note at the bottom of this function for what to do
#' if the boundary is curved/irregular instead.
#'
#' REQUIRES an interactive R session with a graphics device that supports
#' mouse clicks (RStudio's Plots pane, or a windows()/X11() device) - this
#' will NOT work when run non-interactively (e.g. via Rscript on a server
#' with no display).
#'
#' @param obj Seurat spatial object containing both tissue pieces
#' @param tiff_path path to the raw full-resolution TIFF (same one used in
#'   plot_split_boundary_on_tiff())
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @param max_dim downsample the TIFF for display/clicking speed (doesn't
#'   affect the accuracy of the final split - clicked points are converted
#'   back to full-res coordinates automatically)
#' @return the input obj, with a "tissue_piece" metadata column set (1 or 2;
#'   there is no "noise" category with this method - every spot gets assigned
#'   to whichever side of the line it falls on)
manual_split_via_click <- function(obj, tiff_path, image_name = NULL, max_dim = 2000) {
  
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("Package 'tiff' is required. Install it with: install.packages('tiff')")
  }
  if (!file.exists(tiff_path)) {
    stop("TIFF not found: ", tiff_path)
  }
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  cat("Reading TIFF...\n")
  img <- tiff::readTIFF(tiff_path)
  img_h <- dim(img)[1]
  img_w <- dim(img)[2]
  
  coords$plot_row <- coords$imagerow
  coords$plot_col <- coords$imagecol
  scale_factor <- 1
  if (is.finite(max_dim) && max(img_h, img_w) > max_dim) {
    scale_factor <- max_dim / max(img_h, img_w)
    new_h <- max(1, round(img_h * scale_factor))
    new_w <- max(1, round(img_w * scale_factor))
    row_idx <- round(seq(1, img_h, length.out = new_h))
    col_idx <- round(seq(1, img_w, length.out = new_w))
    img <- img[row_idx, col_idx, , drop = FALSE]
    coords$plot_row <- coords$imagerow * scale_factor
    coords$plot_col <- coords$imagecol * scale_factor
    img_h <- new_h
    img_w <- new_w
  }
  
  op <- graphics::par(mar = c(0, 0, 2, 0), mfrow = c(1, 1))
  on.exit(graphics::par(op))
  plot(0, 0, type = "n", xlim = c(0, img_w), ylim = c(img_h, 0), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = "Click 2 points to draw the dividing line between the two tissues")
  graphics::rasterImage(img, 0, img_h, img_w, 0)
  graphics::points(coords$plot_col, coords$plot_row, col = "cyan", pch = 16, cex = 0.3)
  
  cat("\n>>> Click TWO points on the plot window to define a straight dividing line.\n")
  cat(">>> The line should cross through the gap between the two tissues.\n")
  cat(">>> Waiting for clicks...\n")
  pts <- graphics::locator(2)
  if (length(pts$x) < 2) {
    stop("Didn't receive 2 clicks - try again. (Make sure the plot window has focus.)")
  }
  
  # convert clicked (possibly downsampled-for-display) coordinates back to full-res
  x1 <- pts$x[1] / scale_factor; y1 <- pts$y[1] / scale_factor
  x2 <- pts$x[2] / scale_factor; y2 <- pts$y[2] / scale_factor
  
  # side-of-line test (2D cross product sign) - classic point-vs-line classification
  side <- with(coords, (imagecol - x1) * (y2 - y1) - (imagerow - y1) * (x2 - x1))
  
  cluster_col <- setNames(rep(1L, ncol(obj)), colnames(obj))
  cluster_col[coords$barcode] <- ifelse(side > 0, 1L, 2L)
  obj$tissue_piece <- factor(cluster_col)
  
  cat("\nPiece sizes:\n")
  print(table(obj$tissue_piece))
  
  # redraw with the split colored, so you can confirm it looks right immediately
  grDevices::dev.new()
  graphics::par(mar = c(0, 0, 2, 0), mfrow = c(1, 1))
  plot(0, 0, type = "n", xlim = c(0, img_w), ylim = c(img_h, 0), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = "Resulting split - check this before exporting!")
  graphics::rasterImage(img, 0, img_h, img_w, 0)
  split_colors <- c("1" = "#E41A1C", "2" = "#377EB8")
  point_colors <- split_colors[as.character(cluster_col[coords$barcode])]
  graphics::points(coords$plot_col, coords$plot_row, col = point_colors, pch = 16, cex = 0.4)
  graphics::lines(pts$x, pts$y, col = "black", lwd = 2, lty = 2)
  graphics::legend("topright", legend = c("piece 1", "piece 2"),
                   col = split_colors, pch = 16, bty = "o", bg = "white")
  
  cat("\nIf this looks wrong, just re-run manual_split_via_click() and click again.\n")
  cat("If it looks right, export with:\n")
  cat('  split_and_export_dual_tissue(obj, base_sample_name = "...", output_root = "...",\n')
  cat('                                piece_names = c("1" = "...", "2" = "..."))\n')
  
  obj
}

#' MANUAL split by typing in a coordinate cutoff, instead of clicking.
#' Splits every spot based on which side of a line it falls on. The line
#' is defined as: col = intercept + slope * row (axis="col", for left/right
#' splits) or row = intercept + slope * col (axis="row", for top/bottom
#' splits) - slope=0 (the default) gives a perfectly straight cutoff; a
#' nonzero slope lets the line tilt slightly to follow an angled gap
#' between the two tissues.
#'
#' Use plot_manual_cutoff_preview() (below) to see the line overlaid on the
#' image BEFORE committing, and adjust intercept/slope until it sits in the gap.
#'
#' @param obj Seurat spatial object containing both tissue pieces
#' @param cutoff the intercept of the dividing line (full-res pixel space,
#'   same units as pxl_row_in_fullres/pxl_col_in_fullres) - the col (or row)
#'   value of the line at row (or col) = 0
#' @param slope how much the line tilts - e.g. for axis="col", slope=0.1
#'   means the line's col value increases by 0.1 pixels for every 1 pixel
#'   increase in row. Default 0 = perfectly vertical/horizontal.
#' @param axis "col" (splits left/right - most common for side-by-side
#'   tissues) or "row" (splits top/bottom)
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @return the input obj, with a "tissue_piece" metadata column set (1 or 2)
label_tissue_pieces_manual_cutoff <- function(obj, cutoff, slope = 0, axis = c("col", "row"), image_name = NULL) {
  axis <- match.arg(axis)
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  if (axis == "col") {
    threshold <- cutoff + slope * coords$imagerow
    piece <- ifelse(coords$imagecol < threshold, 1L, 2L)
  } else {
    threshold <- cutoff + slope * coords$imagecol
    piece <- ifelse(coords$imagerow < threshold, 1L, 2L)
  }
  
  cluster_col <- setNames(rep(1L, ncol(obj)), colnames(obj))
  cluster_col[coords$barcode] <- piece
  obj$tissue_piece <- factor(cluster_col)
  
  cat("Piece sizes (axis =", axis, ", cutoff =", cutoff, ", slope =", slope, "):\n")
  print(table(obj$tissue_piece))
  
  obj
}

#' Preview a manual coordinate cutoff overlaid on the raw TIFF, before
#' committing to it. Draws a straight (optionally tilted) line, colored
#' dots on either side, so you can adjust cutoff/slope and re-run until
#' the line sits cleanly in the gap between tissues.
#'
#' @param obj Seurat spatial object
#' @param tiff_path path to the raw full-resolution TIFF
#' @param cutoff the line's intercept (same units/meaning as in
#'   label_tissue_pieces_manual_cutoff())
#' @param slope how much the line tilts (0 = straight; see
#'   label_tissue_pieces_manual_cutoff() for the exact meaning) - must match
#'   what you'll pass to label_tissue_pieces_manual_cutoff()
#' @param axis "col" or "row" - must match what you'll pass to
#'   label_tissue_pieces_manual_cutoff()
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @param max_dim downsample the TIFF for display speed (doesn't affect the
#'   actual cutoff/slope values' accuracy)
plot_manual_cutoff_preview <- function(obj, tiff_path, cutoff, slope = 0, axis = c("col", "row"),
                                       image_name = NULL, max_dim = 2000) {
  axis <- match.arg(axis)
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("Package 'tiff' is required. Install it with: install.packages('tiff')")
  }
  if (!file.exists(tiff_path)) stop("TIFF not found: ", tiff_path)
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  img <- tiff::readTIFF(tiff_path)
  img_h_full <- dim(img)[1]; img_w_full <- dim(img)[2]
  img_h <- img_h_full; img_w <- img_w_full
  
  coords$plot_row <- coords$imagerow
  coords$plot_col <- coords$imagecol
  scale_factor <- 1
  if (is.finite(max_dim) && max(img_h, img_w) > max_dim) {
    scale_factor <- max_dim / max(img_h, img_w)
    new_h <- max(1, round(img_h * scale_factor)); new_w <- max(1, round(img_w * scale_factor))
    row_idx <- round(seq(1, img_h, length.out = new_h))
    col_idx <- round(seq(1, img_w, length.out = new_w))
    img <- img[row_idx, col_idx, , drop = FALSE]
    coords$plot_row <- coords$imagerow * scale_factor
    coords$plot_col <- coords$imagecol * scale_factor
    img_h <- new_h; img_w <- new_w
  }
  
  if (axis == "col") {
    threshold <- cutoff + slope * coords$imagerow
    piece <- ifelse(coords$imagecol < threshold, "1", "2")
    # line endpoints in full-res space: (col, row) at row=0 and row=img_h_full
    line_row_full <- c(0, img_h_full)
    line_col_full <- cutoff + slope * line_row_full
  } else {
    threshold <- cutoff + slope * coords$imagecol
    piece <- ifelse(coords$imagerow < threshold, "1", "2")
    # line endpoints in full-res space: (col, row) at col=0 and col=img_w_full
    line_col_full <- c(0, img_w_full)
    line_row_full <- cutoff + slope * line_col_full
  }
  split_colors <- c("1" = "#E41A1C", "2" = "#377EB8")
  
  op <- graphics::par(mar = c(0, 0, 2, 0), mfrow = c(1, 1))
  on.exit(graphics::par(op))
  title_txt <- paste0("Preview: ", axis, " = ", cutoff,
                      if (slope != 0) paste0(" + ", slope, " * ", if (axis == "col") "row" else "col") else "")
  plot(0, 0, type = "n", xlim = c(0, img_w), ylim = c(img_h, 0), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = title_txt)
  graphics::rasterImage(img, 0, img_h, img_w, 0)
  graphics::points(coords$plot_col, coords$plot_row, col = split_colors[piece], pch = 16, cex = 0.4)
  graphics::lines(line_col_full * scale_factor, line_row_full * scale_factor,
                  col = "black", lwd = 2, lty = 2)
  graphics::legend("topright", legend = c("piece 1", "piece 2"),
                   col = split_colors, pch = 16, bty = "o", bg = "white")
  
  cat("Piece sizes at this cutoff:\n")
  print(table(piece))
  cat("If the dashed line isn't in the gap yet, try a different 'cutoff' or 'slope' value.\n")
}

#' MANUAL split using MULTIPLE line constraints combined, for shapes a single
#' straight line can't capture - e.g. an L-shaped tissue, or a trapezoidal
#' slice that needs both a column cutoff AND a row cutoff to carve out
#' correctly. Each constraint is one line (same col/row/slope idea as
#' label_tissue_pieces_manual_cutoff()); a spot is "piece 1" only if it
#' satisfies ALL constraints (logical AND) - this naturally produces
#' rectangular, trapezoidal, or L-cornered regions depending on how many
#' constraints you combine and which direction each one points.
#'
#' @param obj Seurat spatial object containing both tissue pieces
#' @param constraints a list of constraints, each a list with:
#'   - axis: "col" or "row"
#'   - cutoff: the line's intercept
#'   - slope: how much the line tilts (default 0)
#'   - direction: "<" (piece 1 = below/left of the line) or ">=" (piece 1 =
#'     above/right of the line) - default "<"
#'   Example for an L-shaped corner cut (piece 1 = everything left of col
#'   2300 AND above row 1800):
#'     list(
#'       list(axis = "col", cutoff = 2300, direction = "<"),
#'       list(axis = "row", cutoff = 1800, direction = "<")
#'     )
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @return the input obj, with a "tissue_piece" metadata column set (1 or 2)
label_tissue_pieces_manual_polygon <- function(obj, constraints, image_name = NULL) {
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  satisfies_all <- rep(TRUE, nrow(coords))
  for (con in constraints) {
    slope <- if (is.null(con$slope)) 0 else con$slope
    direction <- if (is.null(con$direction)) "<" else con$direction
    
    if (con$axis == "col") {
      threshold <- con$cutoff + slope * coords$imagerow
      value <- coords$imagecol
    } else if (con$axis == "row") {
      threshold <- con$cutoff + slope * coords$imagecol
      value <- coords$imagerow
    } else {
      stop("Each constraint's axis must be 'col' or 'row', got: ", con$axis)
    }
    
    ok <- if (direction == "<") value < threshold else value >= threshold
    satisfies_all <- satisfies_all & ok
  }
  
  cluster_col <- setNames(rep(1L, ncol(obj)), colnames(obj))
  cluster_col[coords$barcode] <- ifelse(satisfies_all, 1L, 2L)
  obj$tissue_piece <- factor(cluster_col)
  
  cat("Piece sizes (", length(constraints), " constraint(s) combined):\n", sep = "")
  print(table(obj$tissue_piece))
  
  obj
}

#' Preview multiple line constraints overlaid on the raw TIFF, before
#' committing. Draws each constraint as its own dashed line across the full
#' image, colors spots by the FINAL combined classification (all constraints
#' AND-ed together), so you can adjust any individual line and re-run until
#' the whole shape looks right.
#'
#' @param obj Seurat spatial object
#' @param tiff_path path to the raw full-resolution TIFF
#' @param constraints same format as label_tissue_pieces_manual_polygon()
#' @param image_name which image slot's coordinates to use (defaults to first)
#' @param max_dim downsample the TIFF for display speed
plot_manual_polygon_preview <- function(obj, tiff_path, constraints, image_name = NULL, max_dim = 2000) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("Package 'tiff' is required. Install it with: install.packages('tiff')")
  }
  if (!file.exists(tiff_path)) stop("TIFF not found: ", tiff_path)
  if (is.null(image_name)) image_name <- Images(obj)[1]
  
  coords <- get_full_res_coords(obj, image_name)
  coords <- coords[coords$barcode %in% colnames(obj), , drop = FALSE]
  
  img <- tiff::readTIFF(tiff_path)
  img_h_full <- dim(img)[1]; img_w_full <- dim(img)[2]
  img_h <- img_h_full; img_w <- img_w_full
  
  coords$plot_row <- coords$imagerow
  coords$plot_col <- coords$imagecol
  scale_factor <- 1
  if (is.finite(max_dim) && max(img_h, img_w) > max_dim) {
    scale_factor <- max_dim / max(img_h, img_w)
    new_h <- max(1, round(img_h * scale_factor)); new_w <- max(1, round(img_w * scale_factor))
    row_idx <- round(seq(1, img_h, length.out = new_h))
    col_idx <- round(seq(1, img_w, length.out = new_w))
    img <- img[row_idx, col_idx, , drop = FALSE]
    coords$plot_row <- coords$imagerow * scale_factor
    coords$plot_col <- coords$imagecol * scale_factor
    img_h <- new_h; img_w <- new_w
  }
  
  satisfies_all <- rep(TRUE, nrow(coords))
  line_endpoints <- list()  # for drawing each constraint's line
  for (i in seq_along(constraints)) {
    con <- constraints[[i]]
    slope <- if (is.null(con$slope)) 0 else con$slope
    direction <- if (is.null(con$direction)) "<" else con$direction
    
    if (con$axis == "col") {
      threshold <- con$cutoff + slope * coords$imagerow
      value <- coords$imagecol
      line_row_full <- c(0, img_h_full)
      line_col_full <- con$cutoff + slope * line_row_full
    } else {
      threshold <- con$cutoff + slope * coords$imagecol
      value <- coords$imagerow
      line_col_full <- c(0, img_w_full)
      line_row_full <- con$cutoff + slope * line_col_full
    }
    ok <- if (direction == "<") value < threshold else value >= threshold
    satisfies_all <- satisfies_all & ok
    line_endpoints[[i]] <- list(col = line_col_full, row = line_row_full)
  }
  piece <- ifelse(satisfies_all, "1", "2")
  split_colors <- c("1" = "#E41A1C", "2" = "#377EB8")
  
  op <- graphics::par(mar = c(0, 0, 2, 0), mfrow = c(1, 1))
  on.exit(graphics::par(op))
  plot(0, 0, type = "n", xlim = c(0, img_w), ylim = c(img_h, 0), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = paste0("Preview: ", length(constraints), " constraint(s) combined"))
  graphics::rasterImage(img, 0, img_h, img_w, 0)
  graphics::points(coords$plot_col, coords$plot_row, col = split_colors[piece], pch = 16, cex = 0.4)
  line_cols <- grDevices::rainbow(length(constraints))
  for (i in seq_along(line_endpoints)) {
    graphics::lines(line_endpoints[[i]]$col * scale_factor, line_endpoints[[i]]$row * scale_factor,
                    col = line_cols[i], lwd = 2, lty = 2)
  }
  graphics::legend("topright", legend = c("piece 1", "piece 2", paste0("constraint ", seq_along(constraints))),
                   col = c(split_colors, line_cols), pch = c(16, 16, rep(NA, length(constraints))),
                   lty = c(NA, NA, rep(2, length(constraints))), lwd = c(NA, NA, rep(2, length(constraints))),
                   bty = "o", bg = "white", cex = 0.75)
  
  cat("Piece sizes at this combination:\n")
  print(table(piece))
  cat("Adjust any individual constraint's cutoff/slope/direction and re-run until the shape looks right.\n")
}


# NOTE: if the boundary between your two tissues is genuinely curved/irregular
# (not well-approximated by one straight line), a straight-line click won't
# capture it well. In that case, Seurat's own InteractiveSpatialPlot() function
# (if available in your installed Seurat version) supports freehand lasso
# selection instead of a single line - check ?InteractiveSpatialPlot. That
# returns a vector of selected cell/barcode names directly, which you can use like:
#   piece1_cells <- InteractiveSpatialPlot(obj, image = Images(obj)[1])
#   obj$tissue_piece <- ifelse(colnames(obj) %in% piece1_cells, 1, 2)

#' Full pipeline: label -> subset -> export each piece as its own
#' Space-Ranger-style folder.
#'
#' @param obj Seurat object (already labeled via label_tissue_pieces(),
#'   OR pass an unlabeled object + this will call label_tissue_pieces()
#'   internally using eps/min_pts/min_component_size)
#' @param base_sample_name e.g. "0223" - used to build output names like "0223_piece1"
#' @param output_root where the two "<name>/" folders will be written
#'   (point this at the same folder your Python script's INPUT_FOLDER reads)
#' @param piece_names optional named vector to give pieces meaningful names,
#'   e.g. c("1" = "0223_pre", "2" = "0223_post") once you've visually confirmed
#'   which cluster is which by comparing against the H&E image / marker genes.
#'   If NULL, defaults to "<base_sample_name>_piece<cluster_number>".
#' @param ... passed to label_tissue_pieces() if obj isn't already labeled
split_and_export_dual_tissue <- function(obj, base_sample_name, output_root,
                                         piece_names = NULL, ...) {
  
  if (!"tissue_piece" %in% colnames(obj@meta.data)) {
    cat("No 'tissue_piece' column found - running label_tissue_pieces() with defaults.\n")
    cat("(Recommended: run label_tissue_pieces() + plot_tissue_pieces() yourself first\n",
        "to confirm the split looks right before exporting.)\n", sep = "")
    res <- label_tissue_pieces(obj, ...)
    obj <- res$obj
  }
  
  sizes <- table(obj$tissue_piece)
  real_clusters <- names(sizes)[names(sizes) != "0"]
  
  if (length(real_clusters) < 2) {
    stop("Found fewer than 2 real tissue clusters (excluding noise). ",
         "Check plot_tissue_pieces() output and try adjusting 'eps' in label_tissue_pieces().")
  }
  if (length(real_clusters) > 2) {
    warning("Found ", length(real_clusters), " clusters, not 2. Exporting all of them - ",
            "review plot_tissue_pieces() to check whether some should be merged ",
            "(eps too small) or excluded as debris (min_component_size too low).")
  }
  
  n_noise <- sum(obj$tissue_piece == "0")
  if (n_noise > 0) {
    cat(n_noise, "spot(s) fell into the 'noise' cluster (0) and will NOT be included ",
        "in either exported piece.\n", sep = "")
  }
  
  if (is.null(piece_names)) {
    piece_names <- setNames(paste0(base_sample_name, "_piece", real_clusters), real_clusters)
  }
  missing_names <- setdiff(real_clusters, names(piece_names))
  if (length(missing_names) > 0) {
    stop("piece_names is missing an entry for cluster(s): ", paste(missing_names, collapse = ", "))
  }
  
  pieces <- list()
  for (cl in real_clusters) {
    sample_name <- piece_names[[cl]]
    cells_in_cluster <- colnames(obj)[obj$tissue_piece == cl]
    cat("\n--- Subsetting cluster", cl, "->", sample_name,
        "(", length(cells_in_cluster), "spots ) ---\n")
    piece_obj <- subset(obj, cells = cells_in_cluster)
    write_seurat_to_spaceranger(piece_obj, sample_name = sample_name, output_root = output_root)
    pieces[[sample_name]] <- piece_obj
  }
  
  cat("\nSplit complete. Exported", length(pieces), "independent sample folder(s) to:", output_root, "\n")
  invisible(pieces)
}

# ============================================================
# Example usage:
#
#   source("split_dual_tissue.R")
#
#   obj <- load_spatial_object("path/to/0223_sample_folder_or_rda", "0223")
#
#   # Step 1: label + ALWAYS visually check before exporting
#   res <- label_tissue_pieces(obj)          # try auto eps first
#   plot_tissue_pieces(res$obj)              # dots-on-Seurat-image view
#   # OR, for a clearer boundary line on the actual raw TIFF:
#   plot_split_boundary_on_tiff(res$obj, tiff_path = "images/Frame_C_S21-003414_and_S19-005879.tif")
#   # if the split looks wrong, retry with an explicit eps, e.g.:
#   # res <- label_tissue_pieces(obj, eps = 400)
#   # plot_split_boundary_on_tiff(res$obj, tiff_path = "images/Frame_C_S21-003414_and_S19-005879.tif")
#
#   # STILL not splitting cleanly? Do it manually by clicking the image:
#   obj <- manual_split_via_click(obj, tiff_path = "images/Frame_C_S21-003414_and_S19-005879.tif")
#   # (shows the result automatically in a new plot window - re-run and click
#   # again if the line was off; this SETS obj$tissue_piece directly, so skip
#   # straight to Step 2 below once it looks right)
#
#   # OR type in a coordinate cutoff instead of clicking (e.g. for left/right
#   # tissues, use axis="col"; iterate on the number until the dashed line
#   # sits in the gap). Add a small 'slope' if the gap is angled rather than
#   # perfectly vertical/horizontal:
#   plot_manual_cutoff_preview(obj, tiff_path = "images/Frame_C_S21-003414_and_S19-005879.tif",
#                               cutoff = 2200, slope = 0.05, axis = "col")
#   # once the preview looks right:
#   obj <- label_tissue_pieces_manual_cutoff(obj, cutoff = 2200, slope = 0.05, axis = "col")
#
#   # For an L-shaped tissue that needs a trapezoidal/corner cut - combine
#   # a column cutoff AND a row cutoff (piece 1 = everything satisfying BOTH):
#   my_constraints <- list(
#     list(axis = "col", cutoff = 2300, slope = 0.01, direction = "<"),
#     list(axis = "row", cutoff = 1800, direction = "<")
#   )
#   plot_manual_polygon_preview(obj, tiff_path = "images/Frame_C_S21-003414_and_S19-005879.tif",
#                                constraints = my_constraints)
#   # adjust any constraint's cutoff/slope/direction, re-run the preview, then commit:
#   obj <- label_tissue_pieces_manual_polygon(obj, constraints = my_constraints)
#
#   # Step 2: once the plot looks right, export both pieces
#   split_and_export_dual_tissue(
#     res$obj,
#     base_sample_name = "0223",
#     output_root = "CRC_project/GSE326101",
#     piece_names = c("1" = "0223_pre", "2" = "0223_post")  # optional, once you know which is which
#   )
#
# After this, "CRC_project/GSE326101/0223_pre/" and ".../0223_post/" are
# two completely independent Space-Ranger-style folders. Your Python
# cell2location script just sees them as two ordinary samples - no changes needed.
# ============================================================