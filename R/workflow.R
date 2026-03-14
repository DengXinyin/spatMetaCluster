#' Run a complete spatial metabolomics clustering workflow
#'
#' This high-level function performs a complete workflow for spatial metabolomics
#' clustering:
#' \enumerate{
#'   \item Extract spectra matrix and pixel metadata from a Cardinal MSI object
#'   \item Remove constant features
#'   \item Apply feature scaling using min-max normalization
#'   \item Run Python UMAP through \pkg{reticulate}
#'   \item Perform clustering on UMAP coordinates
#'   \item Build a clustering result table
#'   \item Write cluster labels back into \code{pixelData(msi_obj)}
#' }
#'
#' This workflow requires a Python environment configured through
#' \pkg{reticulate} when UMAP is run via the Python backend.
#'
#' @param msi_obj A Cardinal MSI object.
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#' @param clustering_method Character string specifying the clustering method.
#' Currently only \code{"kmeans"} is supported. Default is \code{"kmeans"}.
#' @param centers Integer; number of clusters when
#' \code{clustering_method = "kmeans"}. Default is \code{10L}.
#' @param metric Character; UMAP distance metric. Must be one of
#' \code{"cosine"} or \code{"euclidean"}. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors for UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of UMAP dimensions. Any positive integer
#' is allowed, but \code{2} or \code{3} is recommended for most visualization
#' and exploratory analysis tasks. Default is \code{2L}.
#' @param umap_seed Integer; random seed for UMAP. Default is \code{2025L}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{spectra_filtered}{The spectra matrix after removing constant features.}
#'   \item{spectra_scaled}{The spectra matrix after feature scaling using
#'   min-max normalization.}
#'   \item{umap_df}{A data.frame containing UMAP coordinates.}
#'   \item{clustering_result}{A list returned by \code{run_clustering()},
#'   containing the raw clustering result and cluster assignments.}
#'   \item{cluster_df}{A data.frame containing embedding, cluster labels,
#'   and metadata.}
#'   \item{msi_obj}{The updated Cardinal MSI object with cluster labels in
#'   pixelData.}
#' }
#'
#' @examples
#' \dontrun{
#' library(Cardinal)
#' msi_obj <- readImzML("example.imzML")
#'
#' res <- spatial_clustering_workflow(
#'   msi_obj = msi_obj,
#'   python_path = "/path/to/python",
#'   clustering_method = "kmeans",
#'   centers = 10L
#' )
#'
#' head(res$umap_df)
#' table(res$cluster_df$cluster)
#' }
#'
#' @export
spatial_clustering_workflow <- function(
    msi_obj,
    python_path = NULL,
    clustering_method = "kmeans",
    centers = 10L,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    umap_seed = 2025L
) {
  extracted <- extract_spectra_matrix(msi_obj)

  spectra <- extracted$spectra
  pixel_info <- extracted$pixel_info

  spectra <- remove_constant_features(spectra)
  spectra_scaled <- apply_feature_scaling(spectra, method = "minmax")

  umap_df <- run_umap_py(
    x = spectra_scaled,
    python_path = python_path,
    metric = metric,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    random_state = umap_seed
  )

  clustering_res <- run_clustering(
    embedding = umap_df,
    method = clustering_method,
    centers = centers
  )

  cluster_df <- build_cluster_dataframe(
    embedding = umap_df,
    cluster = clustering_res$cluster,
    pixel_info = pixel_info
  )

  msi_obj_updated <- attach_cluster_to_pixeldata(
    msi_obj = msi_obj,
    cluster_df = cluster_df
  )

  list(
    spectra_filtered = spectra,
    spectra_scaled = spectra_scaled,
    umap_df = umap_df,
    clustering_result = clustering_res,
    cluster_df = cluster_df,
    msi_obj = msi_obj_updated
  )
}

