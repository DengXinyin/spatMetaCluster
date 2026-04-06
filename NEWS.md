# ARBO 0.5.1

## New features
- Added two UMAP input strategies in `spatial_clustering_workflow()`:
  - `"scaled"` for min-max scaled spectra
  - `"l2_pca"` for row-wise L2-normalized spectra followed by PCA with `irlba`
- Added `pca_n_components` to control the number of retained principal components
  in the `"l2_pca"` workflow.
- Enforced `metric = "euclidean"` when `umap_input_method = "l2_pca"`.

## Improvements
- Improved compatibility when writing clustering results back to Cardinal
  MSI objects.
- Updated `SEMs_screen()` to better handle different output structures returned
  by `Cardinal::colocalized()`.
- Removed unnecessary use of `magrittr` pipe import.

## Package changes
- Added `irlba` as a package dependency for PCA-based preprocessing.
- Updated workflow documentation and examples to reflect the new UMAP input options.


# ARBO 0.5.0

## New features
- Added two UMAP input strategies in `spatial_clustering_workflow()`:
  - `"scaled"` for min-max scaled spectra
  - `"l2_pca"` for row-wise L2-normalized spectra followed by PCA with `irlba`
- Added `pca_n_components` to control the number of retained principal components
  in the `"l2_pca"` workflow.
- Enforced `metric = "euclidean"` when `umap_input_method = "l2_pca"`.

## Package changes
- Added `irlba` as a package dependency for PCA-based preprocessing.
- Updated workflow documentation and examples to reflect the new UMAP input options.


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
