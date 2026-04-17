# ARBO 1.0.1

## Release highlights

- Added Zenodo DOI information for software citation.
- Added image origin and usage notes for the artificially constructed MSI dataset **Miss Teng**.

## DOI

- Project concept DOI: https://doi.org/10.5281/zenodo.19622206
- Version DOI for ARBO v1.0.1: https://doi.org/10.5281/zenodo.19622207


# ARBO 1.0.0

## Release highlights
- First stable 1.0.0 release of `ARBO`.
- Finalized the main spatial clustering workflow for Cardinal MSI objects.
- Completed package documentation, examples, and vignette materials.
- Consolidated support for two UMAP preprocessing strategies:
  - `"scaled"`: min-max scaled spectra
  - `"l2_pca"`: row-wise L2-normalized spectra followed by PCA with `irlba`
- Stabilized integration of clustering results back into Cardinal MSI objects.
- Finalized utilities for spatially enriched metabolite screening and MSI image visualization.

## Improvements
- Improved robustness of clustering result attachment to MSI objects.
- Refined README and workflow documentation for package release.

# ARBO 0.5.1

## Improvements
- Improved compatibility when writing clustering results back to Cardinal MSI objects.
- Updated `SEMs_screen()` to better handle different output structures returned by `Cardinal::colocalized()`.
- Removed unnecessary use of `magrittr` pipe import.
- Refined workflow documentation and examples.

# ARBO 0.5.0

## New features
- Added two UMAP input strategies in `spatial_clustering_workflow()`:
  - `"scaled"` for min-max scaled spectra
  - `"l2_pca"` for row-wise L2-normalized spectra followed by PCA with `irlba`
- Added `pca_n_components` to control the number of retained principal components in the `"l2_pca"` workflow.
- Enforced `metric = "euclidean"` when `umap_input_method = "l2_pca"`.

## Package changes
- Added `irlba` as a package dependency for PCA-based preprocessing.

# ARBO 0.4.1

## New features
- Added parallel execution support for Python UMAP via `n_jobs`.
- `run_umap_py()` now allows `random_state = NULL` for non-deterministic parallel UMAP.
- `spatial_clustering_workflow()` now supports `n_jobs` and `umap_seed = NULL`.

## Improvements
- Improved documentation for the interaction between UMAP random seeds and parallel execution.
- Clarified that `n_jobs` only affects the Python UMAP step, not downstream clustering in R.

# ARBO 0.4.0
- Renamed the package to `ARBO`.
- Updated package metadata and documentation.
- Refined clustering workflows and result integration into Cardinal MSI objects.

# spatMetaCluster 0.3.1
- Updated vignettes.

# spatMetaCluster 0.3.0
- Added `SEMs_screen()` for spatially enriched metabolite screening.
- Added `image2ggplot()` to reconstruct `Cardinal::image()` output with `ggplot2`.
- Added `msi_img_overlay()` for layered MSI image overlay with internal legend handling.
- Updated `cherry_tomato_msi` with additional m/z annotations and regenerated data.

# spatMetaCluster 0.2.0
- Added a unified clustering interface via `run_clustering()`.
- Renamed `spatial_kmeans_workflow()` to `spatial_clustering_workflow()`.
- Replaced `minmax_normalize()` with `apply_feature_scaling(method = "minmax")`.
- Updated UMAP-related documentation and parameter descriptions.
