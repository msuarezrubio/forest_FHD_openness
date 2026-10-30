# Below-canopy openness and foliage height diversity from portable mobile laser scanning

R code to derive two vegetation structure parameters, below-canopy openness and foliage height diversity (FHD), from height-normalised point clouds of circular forest plots. The code accompanies:

Dammert, L., Arav, R., Suarez-Rubio, M. (2026). Three-dimensional habitat characterisation for forest-dwelling insectivorous bats. Forest Ecology and Management. [DOI]

The point clouds were acquired with a GeoSLAM ZEB Horizon portable mobile laser scanner at 39 sampling points in the Vienna Woods Biosphere Reserve, Austria. The third parameter used in the paper, the number of large-diameter trees, was derived with the treeX workflow (Gollob et al., 2020; Tockner et al., 2022) and is not part of this repository.

# Requirements
R (developed with R 4.5.0), lidR, ggplot2 

# Input data

1. Pre-processed point clouds (.laz, one file per plot) in the folder set as plot_path. Each file must be cropped to a 30 m radius around the plot centre, ground-classified (ground = class 2), and height-normalised, so that Z is the height above ground.

In the paper, ground points were classified with cloth simulation filtering (lidR with RCSF; slope smoothing, class threshold 0.3 m, cloth resolution 0.5 m), and the point clouds were normalised with a 0.3 m digital terrain model interpolated by inverse distance weighting (k = 1000, power = 2). These steps follow the lidR book and are not included here.

The plot ID is taken from the file name as the text before the first underscore (e.g. B01_scan.laz → B01).

2. Tree table (CSV, semicolon-separated, decimal comma) with one row per segmented tree and at least these columns:
Plot:	Plot ID, matching the point cloud file names; est.height:	Estimated tree height (m)

The tallest tree in each plot sets the height of the analysis cylinder.
The point clouds are not included in this repository. [They are available from the corresponding author upon reasonable request.]

# Method

Each plot is processed by process_plot(), which calls the functions below.

1. Voxelisation (get_voxels()). Ground points are removed, and the vegetation point cloud is voxelised at 0.2 m, approximating the smallest wingspan among European insectivorous bats. The voxel grid covers the bounding box of the point cloud and stores the number of points per voxel (0 for empty voxels).

2. Occupancy threshold (find_threshold()). A plot-specific threshold separates occupied from unoccupied voxels to account for differences in point density (following Vonderach et al., 2012). For candidate thresholds of 1–20 points, the number of voxels containing at least that many points is counted. The pair of successive thresholds with the largest decrease in voxel count is identified, and the final threshold is set one point above the higher threshold of that pair.

3. Below-canopy openness. A cylinder with a 30 m radius, extending from the ground to the height of the tallest tree in the plot.

4. Adaptive vertical layering (make_layers()). The first layer comprises all voxels from 0 to 2 m, which corresponds to the shrub layer and to the height of the bat detectors. Above 2 m, the profile is divided into 0.5 m slices. A layer is closed as soon as the proportion of occupied voxels among all voxels of the grid in the slice reaches 0.05; otherwise, the next slice is added until it does. Because the grid covers the square bounding box of the circular plot, this proportion refers to the bounding box. The topmost layer is closed at $h_{max}$, whether or not it reaches the threshold. The threshold of 0.05 was calibrated on the structurally simplest plot, which it resolves into exactly one layer above the shrub layer.

5. Foliage height diversity (FHD) is the Simpson index of the number of occupied voxels per layer.

# Usage
1. In section 1 of the script (Settings), set plot_path, tree_file and out_path.
2. On the first run, keep voxelise <- TRUE. The voxel grid of each plot is then saved as <file>_voxelisation_without_ground.csv next to the .laz file. Later runs can use voxelise <- FALSE, which reads these files and is much faster.
3. Run the whole script. It processes all .laz files in plot_path.


# Citation

If you use this code, please cite the article above and this repository.

