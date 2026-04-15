# ARBO: Rapid spatial clustering and enriched metabolite discovery for MSI data

`ARBO` provides utilities for preprocessing, visualization, and spatial
clustering of metabolomics mass spectrometry imaging (MSI) data, with support
for Python-based UMAP embedding through `reticulate`, multiple clustering
methods, image reconstruction utilities, and integration of clustering results
into Cardinal MSI objects.

The package currently supports two UMAP input strategies in the main workflow:

- **`"scaled"`**: remove constant features, apply min-max scaling, and run UMAP
- **`"l2_pca"`**: remove constant features, perform row-wise L2 normalization,
  reduce dimensionality with PCA using `irlba`, and run UMAP on PCA scores

This allows users to choose either a direct feature-scaled workflow or an
L2-normalized PCA-based workflow before clustering.

## Installation

If installation fails due to missing `Cardinal`, install `Cardinal` first
using Bioconductor, then install `ARBO` from GitHub.

```r
# install.packages("remotes")
remotes::install_github("DengXinyin/ARBO")
```



## R dependencies
The package depends on:
- `Cardinal`
- `reticulate`
- `irlba`

Other imported packages are installed automatically with ARBO.

If Cardinal is not already installed, you may install it with:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("Cardinal")
```

Alternatively, the development version can be installed from GitHub:
```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("kuwisdelu/Cardinal", ref = remotes::github_release())
```

Optional clustering backends used by some methods may require additional
packages such as:
- `mclust`
- `e1071`
- `dbscan`

## Python dependencies
Some functions in ARBO require a Python environment configured for
reticulate, especially for UMAP embedding.

A tested conda environment is:
```bash
conda create -n dxy_python9 python=3.9 -y
conda activate dxy_python9
conda install -c conda-forge numpy=1.24.4 umap-learn=0.5.7 -y
```

You can verify the installed versions with:
```bash
python -c "import numpy; print(numpy.__version__)"
python -c "import umap; print(umap.__version__)"
```

## Notes on compatibility
Other Python versions may also work, but python, numpy, and umap-learn
should be version-compatible. If a newer Python environment is used, please
ensure that the required packages can be successfully imported through
reticulate.

In R, you may point reticulate to the desired Python environment, for example:
```r
reticulate::use_condaenv("dxy_python9", required = TRUE)
```

You can also check whether the Python UMAP backend is available with:
```r
check_umap_python_env()
```

## Example data
The package includes a toy spatial metabolomics imaging dataset:
```r
library(ARBO)
data(cherry_tomato_msi)
cherry_tomato_msi
```

## Example workflow
```r
# 1. Default scaled workflow
library(ARBO)
library(ggplot2)
data(cherry_tomato_msi)

cherry_tomato_msi <- cherry_tomato_msi |>
  Cardinal::peakPick(SNR = 2) |>
  Cardinal::peakAlign()

res <- spatial_clustering_workflow(
  msi_obj = cherry_tomato_msi,
  python_path = "/path/to/conda/env/bin/python",
  clustering_method = "kmeans",
  centers = 2L,
  umap_input_method = "scaled",
  metric = "cosine",
  n_neighbors = 15L,
  min_dist = 0.1,
  n_components = 2L,
  n_jobs = 1L,
  umap_seed = NULL,
  verbose = TRUE
)

cluster_df <- res$cluster_df

ggplot(cluster_df, aes(x, y, color = factor(cluster))) +
  geom_point(size = 1) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  labs(color = "Cluster")
  
# 2. L2-normalized PCA workflow
res_l2_pca <- spatial_clustering_workflow(
  msi_obj = cherry_tomato_msi,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  clustering_method = "kmeans",
  centers = 2L,
  umap_input_method = "l2_pca",
  pca_n_components = 30L,
  metric = "euclidean",
  n_neighbors = 15L,
  min_dist = 0.1,
  n_components = 2L,
  n_jobs = 4L,
  umap_seed = NULL,
  verbose = TRUE
)
  
```
When umap_input_method = "l2_pca", the UMAP metric must be set to
"euclidean".

If you want UMAP parallelism through Python umap-learn, it is recommended to
use umap_seed = NULL together with n_jobs > 1.


## Vignette
A worked example using the cherry tomato toy dataset is available in the package
vignette:
```r
browseVignettes("ARBO")
```

## Project status
`ARBO` is a stable research utility package for spatial metabolomics MSI workflows.
Version 1.0.0 provides a documented and tested workflow for preprocessing,
Python-based UMAP embedding, clustering, visualization, and enriched metabolite screening.

## Main functions
- `extract_spectra_matrix()` – extract spectra and pixel metadata from a Cardinal MSI object
- `remove_constant_features()` – remove zero-variance features
- `apply_feature_scaling()` – scale features
- `run_umap_py()` – run Python UMAP via `reticulate`
- `run_clustering()` – cluster embedding coordinates, e.g. k-means
- `spatial_clustering_workflow()` – complete end-to-end workflow with `"scaled"` and `"l2_pca"` UMAP input strategies
- `add_clusters_to_msi()` – write clustering results back into a Cardinal MSI object
- `SEMs_screen()` – screen spatially enriched metabolites using SSC and colocalization
- `image2ggplot()` – reconstruct `Cardinal::image()` output with `ggplot2`
- `msi_img_overlay()` – overlay multiple MSI images with internal legend handling


## Relationship to Cardinal
`ARBO` is an extension package for workflows based on `Cardinal` MSI objects.
It uses the Cardinal ecosystem for MSI data representation and visualization,
and adds utilities for preprocessing, Python-based UMAP embedding, multiple
clustering methods, and cluster label integration.

Some workflow demonstrations were inspired by the `Cardinal` tutorials.
However, the package code and the included `cherry_tomato_msi` example dataset
were independently created for `ARBO`.
