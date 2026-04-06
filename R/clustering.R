#' Run clustering on embedding coordinates
#'
#' This function performs clustering on low-dimensional embedding coordinates,
#' such as UMAP results. Supported methods include k-means clustering,
#' Gaussian mixture modeling (GMM), fuzzy c-means clustering (FCM),
#' density-based spatial clustering of applications with noise (DBSCAN),
#' and hierarchical clustering.
#'
#' @param embedding A data.frame or matrix containing embedding coordinates.
#' Typically the output of \code{run_umap_py()} or another dimensionality
#' reduction function. All columns must be numeric.
#' @param method Character string specifying the clustering method.
#' Supported options are \code{"kmeans"}, \code{"gmm"}, \code{"fcm"},
#' \code{"dbscan"}, and \code{"hclust"}.
#' @param centers Integer; number of clusters. Used when
#' \code{method = "kmeans"}, \code{"gmm"}, \code{"fcm"}, or \code{"hclust"}.
#' Default is \code{10L}.
#' @param eps Numeric; neighborhood radius used for
#' \code{method = "dbscan"}. Must be provided when using DBSCAN.
#' @param minPts Integer; minimum number of points required to form a dense
#' region for \code{method = "dbscan"}. Default is \code{5L}.
#' @param hclust_method Character string specifying the agglomeration method
#' for hierarchical clustering. Passed to \code{stats::hclust()}.
#' Default is \code{"ward.D2"}.
#' @param seed Integer; random seed used for methods with stochastic behavior,
#' including \code{"kmeans"} and \code{"fcm"}. Default is \code{2026L}.
#' @param nstart Integer; number of random sets used by
#' \code{stats::kmeans()}. Default is \code{10L}.
#'
#' @return A list containing clustering results with at least the following
#' components:
#' \describe{
#'   \item{method}{The clustering method used.}
#'   \item{result}{The raw clustering result object returned by the underlying method.}
#'   \item{cluster}{A factor vector of cluster assignments.}
#' }
#'
#' Additional components may be included depending on the method:
#' \describe{
#'   \item{membership}{For \code{method = "fcm"}, the soft membership matrix.}
#' }
#'
#' For \code{method = "dbscan"}, cluster label \code{0} indicates noise points.
#'
#' @export
run_clustering <- function(
    embedding,
    method = c("kmeans", "gmm", "fcm", "dbscan", "hclust"),
    centers = 10L,
    eps = NULL,
    minPts = 5L,
    hclust_method = "ward.D2",
    seed = 2026L,
    nstart = 10L
) {
  if (!is.data.frame(embedding) && !is.matrix(embedding)) {
    stop("'embedding' must be a data.frame or matrix.")
  }

  if (ncol(embedding) < 2) {
    stop("'embedding' must contain at least 2 columns.")
  }

  embedding_mat <- as.matrix(embedding)

  if (!is.numeric(embedding_mat)) {
    stop("'embedding' must contain only numeric values.")
  }

  if (anyNA(embedding_mat)) {
    stop("'embedding' must not contain missing values.")
  }

  method <- match.arg(method)

  if (!is.numeric(centers) || length(centers) != 1 || is.na(centers)) {
    stop("'centers' must be a single positive integer.")
  }
  centers <- as.integer(centers)
  if (centers < 1L) {
    stop("'centers' must be a positive integer.")
  }

  if (!is.numeric(minPts) || length(minPts) != 1 || is.na(minPts)) {
    stop("'minPts' must be a single positive integer.")
  }
  minPts <- as.integer(minPts)
  if (minPts < 1L) {
    stop("'minPts' must be a positive integer.")
  }

  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("'seed' must be a single integer.")
  }
  seed <- as.integer(seed)

  if (!is.numeric(nstart) || length(nstart) != 1 || is.na(nstart)) {
    stop("'nstart' must be a single positive integer.")
  }
  nstart <- as.integer(nstart)
  if (nstart < 1L) {
    stop("'nstart' must be a positive integer.")
  }

  if (!is.character(hclust_method) || length(hclust_method) != 1 || is.na(hclust_method)) {
    stop("'hclust_method' must be a single character string.")
  }

  if (method %in% c("kmeans", "gmm", "fcm", "hclust")) {
    if (centers > nrow(embedding_mat)) {
      stop("'centers' cannot be greater than the number of rows in 'embedding'.")
    }
  }

  if (method == "kmeans") {
    set.seed(seed)
    km <- stats::kmeans(embedding_mat, centers = centers, nstart = nstart)

    return(list(
      method = method,
      result = km,
      cluster = as.factor(km$cluster)
    ))
  }

  if (method == "gmm") {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("Package 'mclust' is required for method = 'gmm'. Please install it.")
    }

    bic <- mclust::mclustBIC(embedding_mat, G = centers)
    gmm <- mclust::summaryMclustBIC(bic, data = embedding_mat)

    return(list(
      method = method,
      result = gmm,
      cluster = as.factor(gmm$classification)
    ))
  }

  if (method == "fcm") {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package 'e1071' is required for method = 'fcm'. Please install it.")
    }

    set.seed(seed)
    fcm <- e1071::cmeans(embedding_mat, centers = centers)

    return(list(
      method = method,
      result = fcm,
      cluster = as.factor(fcm$cluster),
      membership = fcm$membership
    ))
  }

  if (method == "dbscan") {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required for method = 'dbscan'. Please install it.")
    }

    if (is.null(eps)) {
      stop("'eps' must be provided when method = 'dbscan'.")
    }
    if (!is.numeric(eps) || length(eps) != 1 || is.na(eps) || eps <= 0) {
      stop("'eps' must be a single positive numeric value.")
    }

    db <- dbscan::dbscan(embedding_mat, eps = eps, minPts = minPts)

    return(list(
      method = method,
      result = db,
      cluster = as.factor(db$cluster)
    ))
  }

  if (method == "hclust") {
    if (nrow(embedding_mat) > 5000L) {
      warning(
        "Hierarchical clustering may be slow and memory-intensive for large datasets.",
        call. = FALSE
      )
    }

    d <- stats::dist(embedding_mat)
    hc <- stats::hclust(d, method = hclust_method)
    cl <- stats::cutree(hc, k = centers)

    return(list(
      method = method,
      result = hc,
      cluster = as.factor(cl)
    ))
  }
}


#' Build a clustering result table
#'
#' This function combines low-dimensional embedding coordinates, cluster labels,
#' and pixel metadata into a single data.frame for downstream visualization or
#' export.
#'
#' @param embedding A data.frame or matrix containing embedding coordinates.
#' @param cluster A vector of cluster labels, typically returned by
#'   \code{run_clustering()}.
#' @param pixel_info A data.frame of pixel metadata returned by
#'   \code{extract_spectra_matrix()}.
#'
#' @return A data.frame containing embedding coordinates, cluster assignments,
#' and pixel metadata.
#'
#' @details
#' The number of rows in \code{embedding}, \code{cluster}, and \code{pixel_info}
#' must be identical. The output always contains a column named
#' \code{cluster}.
#'
#' @export
build_cluster_dataframe <- function(embedding, cluster, pixel_info) {
  if (!is.data.frame(embedding) && !is.matrix(embedding)) {
    stop("'embedding' must be a data.frame or matrix.")
  }

  embedding_df <- as.data.frame(embedding)

  if (!is.data.frame(pixel_info)) {
    stop("'pixel_info' must be a data.frame.")
  }

  if (length(cluster) != nrow(embedding_df)) {
    stop("Length of 'cluster' must match the number of rows in 'embedding'.")
  }

  if (nrow(pixel_info) != nrow(embedding_df)) {
    stop("Number of rows in 'pixel_info' must match the number of rows in 'embedding'.")
  }

  if (!"pixel_ID" %in% colnames(pixel_info)) {
    stop("'pixel_info' must contain a 'pixel_ID' column.")
  }

  if (anyNA(pixel_info$pixel_ID)) {
    stop("'pixel_info$pixel_ID' must not contain NA values.")
  }

  if (anyDuplicated(pixel_info$pixel_ID)) {
    stop("'pixel_info$pixel_ID' must be unique.")
  }

  cluster_df <- data.frame(
    pixel_ID = pixel_info$pixel_ID,
    embedding_df,
    cluster = cluster,
    pixel_info[, setdiff(colnames(pixel_info), "pixel_ID"), drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cluster_df
}


#' Attach clustering results back to Cardinal pixelData
#'
#' This function writes cluster labels from a clustering result table back into
#' `pixelData(msi_obj)`.
#'
#' @param msi_obj A Cardinal MSI object.
#' @param cluster_df A data.frame containing at least the columns
#'   `pixel_ID` and `cluster`.
#'
#' @return The input MSI object with updated cluster labels in `pixelData(msi_obj)`.
#' @export
attach_cluster_to_pixeldata <- function(msi_obj, cluster_df) {
  if (!is.data.frame(cluster_df)) {
    stop("'cluster_df' must be a data.frame.")
  }

  required_cols <- c("pixel_ID", "cluster")
  missing_cols <- setdiff(required_cols, colnames(cluster_df))
  if (length(missing_cols) > 0) {
    stop(
      "'cluster_df' must contain the following columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  msi_obj <- ensure_pixel_id(msi_obj)

  pd <- Cardinal::pixelData(msi_obj)
  pd_df <- as.data.frame(pd)

  idx <- match(pd_df$pixel_ID, cluster_df$pixel_ID)

  if (anyNA(idx)) {
    missing_ids <- pd_df$pixel_ID[is.na(idx)]
    stop(
      "Some pixel_ID values in 'pixelData(msi_obj)' were not found in 'cluster_df'. ",
      "Examples: ", paste(utils::head(missing_ids, 5), collapse = ", ")
    )
  }

  pd$cluster <- cluster_df$cluster[idx]
  Cardinal::pixelData(msi_obj) <- pd

  msi_obj
}
