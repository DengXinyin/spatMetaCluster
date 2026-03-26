# ARBO: Rapid spatial clustering and enriched metabolite discovery for MSI data
`ARBO` provides utilities for preprocessing, visualization, and spatial
clustering of metabolomics mass spectrometry imaging (MSI) data, with support
for Python-based UMAP embedding through `reticulate`, multiple clustering
methods, image reconstruction utilities, and integration of clustering results
into Cardinal MSI objects.

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

Additional packages used for vignettes and examples may include:
- `knitr`
- `rmarkdown`
- `ggplot2`

If `Cardinal` is not already installed, you may install it with:
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

## Python dependencies
Some functions in `ARBO` require a Python environment configured for
`reticulate`, especially for UMAP embedding.

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
Other Python versions may also work, but `python`, `numpy`, and `umap-learn`
should be version-compatible. If a newer Python environment is used, please
ensure that the required packages can be successfully imported through
`reticulate`.

In R, you may point `reticulate` to the desired Python environment, for example:
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
library(ARBO)
library(ggplot2)
data(cherry_tomato_msi)

cherry_tomato_msi <- cherry_tomato_msi |>
  Cardinal::peakPick(SNR = 2) |>
  Cardinal::peakAlign()

res <- spatial_clustering_workflow(
  msi_obj = cherry_tomato_msi,
  python_path = "/home/hsinyinteng/miniconda3/envs/dxy_python9/bin/python",
  clustering_method = "kmeans",
  centers = 2L,  
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
  
```

## Vignette
A worked example using the cherry tomato toy dataset is available in the package
vignette:
```r
browseVignettes("ARBO")
```

## Project status
`ARBO` is a lightweight research utility package for spatial
metabolomics clustering workflows. The current version is functional and
documented. Future extensions may include additional preprocessing or clustering
options.

## Main functions
- `extract_spectra_matrix()` – extract spectra and pixel metadata from a Cardinal MSI object
- `remove_constant_features()` – remove zero-variance features
- `apply_feature_scaling()` – scale features
- `run_umap_py()` – run Python UMAP via `reticulate`
- `run_clustering()` – cluster embedding coordinates, e.g. k-means
- `spatial_clustering_workflow()` – complete end-to-end workflow

## Relationship to Cardinal
`ARBO` is an extension package for workflows based on
`Cardinal` MSI objects. It uses the Cardinal ecosystem for MSI data
representation and visualization, and adds utilities for preprocessing,
Python-based UMAP embedding, k-means clustering, and cluster label
integration.

Some workflow demonstrations were inspired by the Cardinal tutorials.
However, the package code and the included `cherry_tomato_msi` example
dataset were independently created for `ARBO`.
