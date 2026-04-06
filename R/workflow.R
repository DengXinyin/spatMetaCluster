# Internal helper: L2-normalize rows of a matrix
.l2_normalize_rows <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  check_spectra_matrix(x)

  row_norms <- sqrt(rowSums(x^2, na.rm = TRUE))
  row_norms[!is.finite(row_norms) | row_norms == 0] <- 1

  x_norm <- x / row_norms
  rownames(x_norm) <- rownames(x)
  colnames(x_norm) <- colnames(x)

  x_norm
}

# Internal helper: run PCA using irlba::prcomp_irlba
.run_pca_irlba <- function(x, n_components = 30L, center = TRUE, scale. = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  check_spectra_matrix(x)

  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required but not installed.")
  }

  if (!is.numeric(n_components) || length(n_components) != 1L || is.na(n_components)) {
    stop("'pca_n_components' must be a single positive integer.")
  }
  n_components <- as.integer(n_components)

  if (n_components < 1L) {
    stop("'pca_n_components' must be >= 1.")
  }

  max_rank <- min(nrow(x) - 1L, ncol(x) - 1L)
  if (max_rank < 1L) {
    stop("Input matrix is too small for PCA after preprocessing.")
  }

  if (n_components > max_rank) {
    stop(
      "'pca_n_components' must be <= min(nrow(x) - 1, ncol(x) - 1). ",
      "Current maximum allowed value is ", max_rank, "."
    )
  }

  pca <- irlba::prcomp_irlba(
    x = x,
    n = n_components,
    center = center,
    scale. = scale.
  )

  scores <- as.data.frame(pca$x)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  rownames(scores) <- rownames(x)

  list(
    model = pca,
    scores = scores
  )
}


#' Run a complete spatial metabolomics clustering workflow
#'
#' This high-level function performs a complete workflow for spatial
#' metabolomics clustering:
#' \enumerate{
#'   \item Extract the spectra matrix and pixel metadata from a Cardinal MSI object
#'   \item Remove constant features
#'   \item Prepare the UMAP input using one of two strategies:
#'   \enumerate{
#'     \item \code{"scaled"}: apply feature scaling using min-max normalization
#'     \item \code{"l2_pca"}: apply row-wise L2 normalization, then PCA using
#'     \code{irlba::prcomp_irlba()}
#'   }
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
#' Two UMAP input strategies are supported:
#' \itemize{
#'   \item \code{"scaled"}: the default general-purpose workflow. After removing
#'   constant features, the spectra matrix is min-max scaled and passed directly
#'   to UMAP.
#'   \item \code{"l2_pca"}: a workflow intended for L2-normalized PCA embedding
#'   followed by UMAP with Euclidean distance. In this mode, the spectra matrix
#'   is row-wise L2-normalized, reduced by PCA using
#'   \code{irlba::prcomp_irlba()}, and the resulting PCA scores are used as input
#'   to UMAP. For this strategy, \code{metric} must be \code{"euclidean"}.
#' }
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
#' @param umap_input_method Character string specifying how the input to UMAP
#' is constructed. Supported options are \code{"scaled"} and \code{"l2_pca"}.
#' Default is \code{"scaled"}.
#'
#' \code{"scaled"} is the default general-purpose workflow: after constant
#' feature removal, min-max scaled spectra are used directly as input to UMAP.
#'
#' \code{"l2_pca"} is intended for L2-normalized PCA embedding followed by
#' Euclidean UMAP: after constant feature removal, spectra are row-wise
#' L2-normalized, reduced by PCA using \code{irlba::prcomp_irlba()}, and the
#' resulting PCA scores are used as input to UMAP. When
#' \code{umap_input_method = "l2_pca"}, \code{metric} must be
#' \code{"euclidean"}.
#' @param pca_n_components Integer; number of principal components retained
#' when \code{umap_input_method = "l2_pca"}. Default is \code{30L}. Ignored
#' when \code{umap_input_method = "scaled"}.
#' @param metric Character; UMAP distance metric. Supported values are
#' \code{"cosine"}, \code{"correlation"}, \code{"euclidean"},
#' \code{"chebyshev"}, \code{"manhattan"}, \code{"minkowski"},
#' \code{"canberra"}, \code{"braycurtis"}, \code{"hamming"}, and
#' \code{"jaccard"}. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors for UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of UMAP dimensions. Default is \code{2L}.
#' @param umap_seed Optional integer; random seed for UMAP. Default is
#' \code{2025L}. Use \code{NULL} to allow non-deterministic execution.
#' @param n_jobs Integer; number of parallel jobs used by the Python UMAP
#' backend. Default is \code{1L}. Use \code{-1L} to request all available CPU
#' cores if supported by the Python backend.
#' @param verbose Logical; whether to print UMAP progress information.
#' Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{spectra_filtered}{A numeric matrix after removing constant features.}
#'   \item{spectra_scaled}{A numeric matrix after feature scaling using
#'   min-max normalization when \code{umap_input_method = "scaled"}; otherwise
#'   \code{NULL}.}
#'   \item{spectra_l2}{A numeric matrix after row-wise L2 normalization when
#'   \code{umap_input_method = "l2_pca"}; otherwise \code{NULL}.}
#'   \item{pca_result}{A list containing the PCA model and PCA scores when
#'   \code{umap_input_method = "l2_pca"}; otherwise \code{NULL}.}
#'   \item{umap_input}{The object used as direct input to UMAP.}
#'   \item{umap_df}{A data.frame containing UMAP coordinates.}
#'   \item{clustering_result}{A list returned by \code{run_clustering()}.}
#'   \item{cluster_df}{A data.frame containing UMAP coordinates, cluster labels,
#'   and pixel metadata.}
#'   \item{msi_obj}{The updated Cardinal MSI object with cluster labels added
#'   to \code{pixelData}.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- spatial_clustering_workflow(
#'   msi_obj = tomato,
#'   python_path = "/path/to/python",
#'   clustering_method = "kmeans",
#'   centers = 2L,
#'   umap_input_method = "scaled",
#'   metric = "cosine",
#'   n_neighbors = 15L,
#'   min_dist = 0.1,
#'   n_components = 2L,
#'   n_jobs = 1L,
#'   umap_seed = NULL,
#'   verbose = TRUE
#' )
#' }
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
    umap_input_method = c("scaled", "l2_pca"),
    pca_n_components = 30L,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    umap_seed = 2025L,
    n_jobs = 1L,
    verbose = TRUE
) {
  clustering_method <- match.arg(clustering_method)
  umap_input_method <- match.arg(umap_input_method)

  if (clustering_method == "dbscan" && is.null(eps)) {
    stop("'eps' must be provided when clustering_method = 'dbscan'.")
  }

  if (umap_input_method == "l2_pca" && metric != "euclidean") {
    stop("When 'umap_input_method = \"l2_pca\"', 'metric' must be 'euclidean'.")
  }

  msi_obj <- ensure_pixel_id(msi_obj)

  extracted <- extract_spectra_matrix(msi_obj)
  msi_obj <- extracted$msi_obj

  spectra <- extracted$spectra
  pixel_info <- extracted$pixel_info

  spectra <- as.matrix(spectra)
  storage.mode(spectra) <- "numeric"
  check_spectra_matrix(spectra)

  spectra <- remove_constant_features(spectra)
  spectra <- as.matrix(spectra)
  storage.mode(spectra) <- "numeric"
  check_spectra_matrix(spectra)

  spectra_scaled <- NULL
  spectra_l2 <- NULL
  pca_result <- NULL
  umap_input <- NULL

  if (umap_input_method == "scaled") {
    spectra_scaled <- apply_feature_scaling(spectra, method = "minmax")
    spectra_scaled <- as.matrix(spectra_scaled)
    storage.mode(spectra_scaled) <- "numeric"
    umap_input <- spectra_scaled
  }

  if (umap_input_method == "l2_pca") {
    spectra_l2 <- .l2_normalize_rows(spectra)
    pca_result <- .run_pca_irlba(
      x = spectra_l2,
      n_components = pca_n_components,
      center = TRUE,
      scale. = FALSE
    )
    umap_input <- as.matrix(pca_result$scores)
    storage.mode(umap_input) <- "numeric"
  }

  umap_df <- run_umap_py(
    x = umap_input,
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
    spectra_l2 = spectra_l2,
    pca_result = pca_result,
    umap_input = umap_input,
    umap_df = umap_df,
    clustering_result = clustering_res,
    cluster_df = cluster_df,
    msi_obj = msi_obj_updated
  )
}
