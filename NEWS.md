# ARBO 0.4.1
## New features
- Added parallel execution support for Python UMAP via `n_jobs`.
- `run_umap_py()` now allows `random_state = NULL` for non-deterministic parallel UMAP.
- `spatial_clustering_workflow()` now supports `n_jobs` and `umap_seed = NULL`.

## Improvements
- Improved documentation for the interaction between UMAP random seeds and parallel execution.
- Clarified that `n_jobs` only affects the Python UMAP step, not downstream clustering in R.

# ARBO 0.4.0
- Rename the package to `ARBO`.
- Update package metadata and documentation.
- Refine clustering workflows and result integration into Cardinal MSI objects.

# spatMetaCluster 0.3.1
- Update vignettes.

# spatMetaCluster 0.3.0
- Add `SEMs_screen()` for spatially enriched metabolite screening.
- Add `image2ggplot()` to reconstruct `Cardinal::image()` output with `ggplot2`.
- Add `msi_img_overlay()` for layered MSI image overlay with internal legend handling.
- Update `cherry_tomato_msi` with additional m/z annotations and regenerated data.

# spatMetaCluster 0.2.0
- Add a unified clustering interface via `run_clustering()`.
- Rename `spatial_kmeans_workflow()` to `spatial_clustering_workflow()`.
- Replace `minmax_normalize()` with `apply_feature_scaling(method = "minmax")`.
- Update UMAP-related documentation and parameter descriptions.
