# =============================================================================
# Below-canopy openness and foliage height diversity (FHD)
# from portable mobile laser scanning point clouds
#
# Accompanies:
#   Dammert, L., Arav, R., Suarez-Rubio, M. (2026). Three-dimensional habitat
#   characterisation for forest-dwelling insectivorous bats.
#   Forest Ecology and Management. [DOI to be added]
#
# Input (see README):
#   - Height-normalised, ground-classified point clouds (.laz), one per plot,
#     cropped to a 30 m radius around the plot centre. The plot ID is the part
#     of the file name before the first underscore (e.g. "B01_scan.laz").
#   - A table of segmented trees with the columns "Plot" and "est.height".
#
# Output (written to `out_path`):
#   - openness_fhd_all_plots.csv: one row per plot
#   - <plot>_threshold_search.csv: occupied voxels per candidate threshold
#   - <plot>_layers.csv: vegetation layers with occupied/unoccupied voxels
#   - <plot>_layer_profile.png: diagnostic vertical profile (optional)
#
# Developed with R 4.5.0, lidR [x.y.z] and ggplot2 [x.y.z].
# =============================================================================

library(lidR)
library(ggplot2)


# ---- 1. Settings ------------------------------------------------------------

plot_path <- "path/to/point_clouds"     # folder containing the .laz files
tree_file <- "path/to/tree_table.csv"   # segmented trees (";"-separated, decimal comma)
out_path  <- "output"                   # folder for the results

voxelise   <- TRUE  # TRUE: voxelise the point clouds and save the voxel grids;
                    # FALSE: read voxel grids saved in a previous run (faster)
make_plots <- TRUE  # write a diagnostic vertical profile for each plot

res             <- 0.2    # voxel size (m)
plot_radius     <- 30     # plot radius (m)
threshold_range <- 1:20   # candidate occupancy thresholds (points per voxel)
shrub_layer_top <- 2      # upper limit of the first vegetation layer (m)
slice_height    <- 0.5    # height of the horizontal slices (m)
layer_threshold <- 0.05   # proportion of occupied voxels needed to close a layer


# ---- 2. Functions -----------------------------------------------------------

# Voxelise the vegetation points of one plot, or read a saved voxel grid.
# Returns all voxels of the grid (the bounding box of the point cloud) with the
# number of points per voxel in column V1 (0 for empty voxels). The grid is
# saved as CSV next to the .laz file and always read back from it, so results
# are identical whether or not the voxelisation is rerun.
get_voxels <- function(laz_file, voxelise) {
  csv_file <- sub("\\.laz$", "_voxelisation_without_ground.csv", laz_file)

  if (voxelise) {
    las <- readLAS(laz_file, select = "xyzc")
    las <- filter_poi(las, Classification != 2L)   # remove ground points
    vox <- voxel_metrics(las, ~length(Z), res = res, all_voxels = TRUE)
    vox[is.na(vox)] <- 0
    write.csv2(vox, csv_file)
  }

  read.csv2(csv_file, fileEncoding = "UTF-8")
}


# Plot-specific occupancy threshold (following Vonderach et al., 2012).
# For each candidate threshold t, the voxels containing at least t points are
# counted. The pair of successive thresholds (t - 1, t) with the largest
# decrease in this count is identified, and the final threshold is t + 1.
find_threshold <- function(points_per_voxel, candidates = threshold_range) {
  n_voxels <- sapply(candidates, function(t) sum(points_per_voxel >= t))
  decrease <- c(NA, abs(diff(n_voxels)))
  t_upper  <- candidates[which.max(decrease)]

  list(threshold = t_upper + 1,
       table = data.frame(threshold = candidates,
                          n_voxels_occupied = n_voxels,
                          decrease = decrease))
}


# Adaptive vertical layering.
# Layer 1 contains all voxels from 0 m up to `shrub_layer_top`. Above, the
# profile is divided into slices of `slice_height`. A layer is closed as soon
# as the proportion of occupied voxels among all voxels of the grid in the
# slice reaches `layer_threshold`; otherwise, the next slice is added. The
# topmost layer is closed at `z_max` (height of the tallest tree), whether or
# not it reaches the threshold. Because the voxel grid spans the square
# bounding box of the circular plot, proportions refer to that box.
make_layers <- function(z_occupied, z_unoccupied, z_max) {
  count_voxels <- function(lower, upper) {
    c(occupied   = sum(z_occupied   >= lower & z_occupied   < upper),
      unoccupied = sum(z_unoccupied >= lower & z_unoccupied < upper))
  }

  bounds <- list(c(0, shrub_layer_top))
  lower  <- shrub_layer_top

  while (lower < z_max) {
    upper <- lower + slice_height
    repeat {
      n <- count_voxels(lower, upper)
      prop_occupied <- if (sum(n) > 0) n[["occupied"]] / sum(n) else 0
      if (prop_occupied >= layer_threshold || upper >= z_max) break
      upper <- upper + slice_height
    }
    bounds[[length(bounds) + 1]] <- c(lower, upper)
    lower <- upper
  }

  layers <- do.call(rbind, lapply(bounds, function(b) {
    n <- count_voxels(b[1], b[2])
    data.frame(lower_m = b[1], upper_m = b[2],
               n_occupied = n[["occupied"]], n_unoccupied = n[["unoccupied"]])
  }))
  n_total <- layers$n_occupied + layers$n_unoccupied
  layers$prop_occupied <- ifelse(n_total > 0, layers$n_occupied / n_total, 0)
  layers
}


# Foliage height diversity: Simpson index of the number of occupied voxels
# per layer, FHD = 1 - sum((n_i / N)^2), with N = sum(n_i).
simpson_fhd <- function(n) {
  p <- n / sum(n)
  1 - sum(p^2)
}


# Diagnostic plot: occupied voxels per vegetation layer over height.
plot_profile <- function(layers, plotid, threshold) {
  p <- ggplot(layers) +
    geom_rect(aes(xmin = 0, xmax = n_occupied, ymin = lower_m, ymax = upper_m),
              fill = "springgreen4", colour = "white") +
    labs(title = paste0("Plot ", plotid,
                        " (occupancy threshold: ", threshold, " points)"),
         x = "Occupied voxels per layer", y = "Height above ground (m)") +
    theme_light()
  ggsave(file.path(out_path, paste0(plotid, "_layer_profile.png")), p,
         width = 6, height = 6, dpi = 300)
}


# Process one plot: occupancy threshold, below-canopy openness, layers and FHD.
process_plot <- function(laz_file) {
  plotid <- strsplit(basename(laz_file), "_")[[1]][1]
  message("Processing plot ", plotid)

  vox <- get_voxels(laz_file, voxelise)

  # Occupied voxels
  th <- find_threshold(vox$V1)
  occupied <- vox$V1 >= th$threshold
  n_occupied <- sum(occupied)

  # Below-canopy openness: unoccupied share of a cylinder with radius
  # `plot_radius`, from the ground to the height of the tallest tree
  heights <- trees$est.height[trees$Plot == plotid]
  if (length(heights) == 0) {
    stop("No trees found for plot ", plotid, " in the tree table.")
  }
  z_max <- max(heights)
  vol_cylinder     <- pi * plot_radius^2 * z_max
  n_voxels_cylinder <- vol_cylinder / res^3
  vol_unoccupied   <- (n_voxels_cylinder - n_occupied) * res^3
  openness         <- vol_unoccupied / vol_cylinder

  # Vertical layering and FHD
  layers <- make_layers(vox$Z[occupied], vox$Z[!occupied], z_max)
  fhd <- simpson_fhd(layers$n_occupied)

  # Plot-level outputs
  write.csv2(th$table, row.names = FALSE,
             file.path(out_path, paste0(plotid, "_threshold_search.csv")))
  write.csv2(layers, row.names = FALSE,
             file.path(out_path, paste0(plotid, "_layers.csv")))
  if (make_plots) plot_profile(layers, plotid, th$threshold)

  data.frame(plot                 = plotid,
             occupancy_threshold  = th$threshold,
             max_tree_height_m    = z_max,
             n_voxels_cylinder    = n_voxels_cylinder,
             n_voxels_occupied    = n_occupied,
             volume_cylinder_m3   = vol_cylinder,
             volume_unoccupied_m3 = vol_unoccupied,
             below_canopy_openness = openness,
             n_layers             = nrow(layers),
             fhd                  = fhd)
}


# ---- 3. Run for all plots ---------------------------------------------------

dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
trees <- read.csv(tree_file, sep = ";", dec = ",")

laz_files <- list.files(plot_path, pattern = "\\.laz$", full.names = TRUE)
results <- do.call(rbind, lapply(laz_files, process_plot))

write.csv2(results, file.path(out_path, "openness_fhd_all_plots.csv"),
           row.names = FALSE)
print(results)
