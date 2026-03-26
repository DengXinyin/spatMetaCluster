#' Run a complete spatial metabolomics clustering workflow
#'
#' This high-level function performs a complete workflow for spatial
#' metabolomics clustering:
#' \enumerate{
#'   \item Extract the spectra matrix and pixel metadata from a Cardinal MSI object
#'   \item Remove constant features
#'   \item Apply feature scaling using min-max normalization
#'   \item Run Python UMAP through \pkg{reticulate}
#'   \item Perform clustering on the UMAP embedding
#'   \item Build a clustering result table
#'   \item Attach cluster labels back to \code{pixelData(msi_obj)}
#' }
#'
#' This workflow is designed for spatial metabolomics data stored in a
#' Cardinal-compatible MSI object. UMAP is computed using the Python
#' \pkg{umap-learn} backend via \pkg{reticulate}.
#'
#' @param msi_obj A Cardinal MSI object.
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#' @param clustering_method Character string specifying the clustering method.
#' Supported options are \code{"kmeans"}, \code{"gmm"}, \code{"fcm"},
#' \code{"dbscan"}, and \code{"hclust"}. Default is \code{"kmeans"}.
#' @param centers Integer; number of clusters when
#' \code{clustering_method = "kmeans"}, \code{"gmm"}, \code{"fcm"},
#' or \code{"hclust"}. Default is \code{10L}.
#' @param eps Numeric; neighborhood radius used when
#' \code{clustering_method = "dbscan"}. Must be provided for DBSCAN.
#' Default is \code{NULL}.
#' @param minPts Integer; minimum number of points required to form a dense
#' region when \code{clustering_method = "dbscan"}. Default is \code{5L}.
#' @param hclust_method Character string specifying the agglomeration method
#' for hierarchical clustering. Passed to \code{stats::hclust()}.
#' Default is \code{"ward.D2"}.
#' @param clustering_seed Integer; random seed for clustering methods with
#' stochastic behavior, including \code{"kmeans"} and \code{"fcm"}.
#' Default is \code{2026L}.
#' @param nstart Integer; number of random starts used by
#' \code{stats::kmeans()}. Default is \code{10L}.
#' @param metric Character; UMAP distance metric. Must be one of
#' \code{"cosine"} or \code{"euclidean"}. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors for UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of UMAP dimensions. Any positive integer
#' is allowed, but \code{2} or \code{3} is recommended for most visualization
#' and exploratory analysis tasks. Default is \code{2L}.
#' @param umap_seed Optional integer; random seed for UMAP. Default is
#' \code{2025L}. Use \code{NULL} to allow non-deterministic execution.
#' In Python \pkg{umap-learn}, a fixed seed may disable or limit parallel
#' execution. If \code{n_jobs != 1L} and \code{umap_seed} is not \code{NULL},
#' the downstream UMAP call may switch the Python-side seed to \code{NULL}
#' to allow parallel execution.
#' @param n_jobs Integer; number of parallel jobs used by the Python UMAP
#' backend. Default is \code{1L}. Use \code{-1L} to request all available CPU
#' cores if supported by the Python backend. This argument affects only the
#' Python UMAP step and does not automatically parallelize downstream clustering
#' performed in R.
#' @param verbose Logical; whether to print UMAP progress information.
#' Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{spectra_filtered}{A numeric matrix after removing constant features.}
#'   \item{spectra_scaled}{A numeric matrix after feature scaling using
#'   min-max normalization.}
#'   \item{umap_df}{A data.frame containing UMAP coordinates.}
#'   \item{clustering_result}{A list returned by \code{run_clustering()},
#'   containing the fitted clustering result and cluster assignments.}
#'   \item{cluster_df}{A data.frame containing UMAP coordinates, cluster labels,
#'   and pixel metadata.}
#'   \item{msi_obj}{The updated Cardinal MSI object with cluster labels added
#'   to \code{pixelData}.}
#' }
#'
#' @details
#' The clustering step is performed on the UMAP embedding generated from the
#' scaled spectra matrix.
#'
#' Method-specific clustering arguments are interpreted as follows:
#' \itemize{
#'   \item \code{"kmeans"} uses \code{centers}, \code{clustering_seed}, and \code{nstart}.
#'   \item \code{"gmm"} uses \code{centers}.
#'   \item \code{"fcm"} uses \code{centers} and \code{clustering_seed}.
#'   \item \code{"dbscan"} uses \code{eps} and \code{minPts}.
#'   \item \code{"hclust"} uses \code{centers} and \code{hclust_method}.
#' }
#'
#' The \code{n_jobs} argument only controls parallel execution in the Python
#' \pkg{umap-learn} backend. It does not parallelize feature extraction,
#' feature scaling, or clustering methods executed in R.
#'
#' For Python \pkg{umap-learn}, fixed random seeds and parallel execution may
#' conflict. In practice, when \code{umap_seed} is fixed, the backend may
#' override \code{n_jobs} to \code{1}. To enable actual parallel execution,
#' set \code{umap_seed = NULL}.
#'
#' For \code{clustering_method = "dbscan"}, cluster label \code{0} indicates
#' noise points.
#'
#' @examples
#' \dontrun{
#' library(Cardinal)
#'
#' msi_obj <- readImzML("example.imzML")
#'
#' ## Reproducible single-thread workflow
#' res1 <- spatial_clustering_workflow(
#'   msi_obj = msi_obj,
#'   python_path = "/path/to/python",
#'   clustering_method = "kmeans",
#'   centers = 10L,
#'   metric = "cosine",
#'   n_neighbors = 10L,
#'   min_dist = 0.05,
#'   n_components = 2L,
#'   umap_seed = 2025L,
#'   n_jobs = 1L
#' )
#'
#' head(res1$umap_df)
#' table(res1$cluster_df$cluster)
#'
#' ## Parallel UMAP workflow
#' ## To allow actual parallelism in Python umap-learn, use umap_seed = NULL
#' res2 <- spatial_clustering_workflow(
#'   msi_obj = msi_obj,
#'   python_path = "/path/to/python",
#'   clustering_method = "kmeans",
#'   centers = 10L,
#'   metric = "cosine",
#'   n_neighbors = 10L,
#'   min_dist = 0.05,
#'   n_components = 2L,
#'   umap_seed = NULL,
#'   n_jobs = 4L
#' )
#'
#' ## DBSCAN example
#' res3 <- spatial_clustering_workflow(
#'   msi_obj = msi_obj,
#'   python_path = "/path/to/python",
#'   clustering_method = "dbscan",
#'   eps = 0.3,
#'   minPts = 10L,
#'   umap_seed = 2025L,
#'   n_jobs = 1L
#' )
#' }
#'
#' @export
spatial_clustering_workflow <- function(
    msi_obj,
    python_path = NULL,
    clustering_method = c("kmeans", "gmm", "fcm", "dbscan", "hclust"),
    centers = 10L,
    eps = NULL,
    minPts = 5L,
    hclust_method = "ward.D2",
    clustering_seed = 2026L,
    nstart = 10L,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    umap_seed = 2025L,
    n_jobs = 1L,
    verbose = TRUE
) {
  clustering_method <- match.arg(clustering_method)

  if (clustering_method == "dbscan" && is.null(eps)) {
    stop("'eps' must be provided when clustering_method = 'dbscan'.")
  }

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
    random_state = umap_seed,
    n_jobs = n_jobs,
    verbose = verbose
  )

  clustering_res <- run_clustering(
    embedding = umap_df,
    method = clustering_method,
    centers = centers,
    eps = eps,
    minPts = minPts,
    hclust_method = hclust_method,
    seed = clustering_seed,
    nstart = nstart
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

