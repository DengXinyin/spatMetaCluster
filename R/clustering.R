#' Run clustering on embedding coordinates
#'
#' This function performs clustering on low-dimensional embedding coordinates,
#' such as UMAP results. Currently, only \code{"kmeans"} clustering is supported,
#' but additional methods may be added in future versions.
#'
#' @param embedding A data.frame or matrix containing embedding coordinates.
#' Typically the output of \code{run_umap_py()}.
#' @param method Character string specifying the clustering method.
#' Currently only \code{"kmeans"} is supported.
#' @param centers Integer; number of clusters. Used when \code{method = "kmeans"}.
#' Default is \code{10L}.
#'
#' @return A list containing clustering results. For \code{method = "kmeans"},
#' the returned list includes:
#' \describe{
#'   \item{method}{The clustering method used.}
#'   \item{result}{The raw \code{stats::kmeans()} result object.}
#'   \item{cluster}{A factor vector of cluster assignments.}
#' }
#'
#' @details
#' For \code{method = "kmeans"}, clustering is performed on all columns of the
#' embedding coordinates, and the random seed is fixed internally at
#' \code{2026} for reproducibility.
#'
#' @examples
#' emb <- data.frame(UMAP1 = rnorm(20), UMAP2 = rnorm(20))
#' res <- run_clustering(emb, method = "kmeans", centers = 3)
#' table(res$cluster)
#'
#' @export
run_clustering <- function(
    embedding,
    method = "kmeans",
    centers = 10L
) {
  if (!is.data.frame(embedding) && !is.matrix(embedding)) {
    stop("'embedding' must be a data.frame or matrix.")
  }

  if (ncol(embedding) < 2) {
    stop("'embedding' must contain at least 2 columns.")
  }

  method <- match.arg(method, choices = c("kmeans"))

  if (!is.numeric(centers) || length(centers) != 1 || is.na(centers)) {
    stop("'centers' must be a single positive integer.")
  }
  centers <- as.integer(centers)
  if (centers < 1L) {
    stop("'centers' must be a positive integer.")
  }

  embedding_mat <- as.matrix(embedding)

  if (!is.numeric(embedding_mat)) {
    stop("'embedding' must contain only numeric values.")
  }

  if (method == "kmeans") {
    set.seed(2026)
    km <- stats::kmeans(embedding_mat, centers = centers)

    return(list(
      method = method,
      result = km,
      cluster = as.factor(km$cluster)
    ))
  }
}


#' Build a clustering result table
#'
#' This function combines low-dimensional embedding coordinates, cluster labels,
#' and pixel metadata into a single data.frame for downstream visualization or
#' export.
#'
#' @param embedding Data frame or matrix with embedding coordinates.
#' @param cluster Cluster labels.
#' @param pixel_info Pixel metadata from \code{extract_spectra_matrix()}.
#'
#' @return A data.frame containing UMAP coordinates, cluster assignments, and
#' selected pixel metadata.
#'
#' @export
build_cluster_dataframe <- function(embedding, cluster, pixel_info) {
  embedding <- as.data.frame(embedding)
  pixel_info <- as.data.frame(pixel_info)

  if (nrow(embedding) != length(cluster)) {
    stop("Length of cluster labels must match number of rows in embedding.")
  }

  if (nrow(pixel_info) != nrow(embedding)) {
    stop("pixel_info must have the same number of rows as embedding.")
  }

  if (!"pixel_ID" %in% colnames(pixel_info)) {
    stop("pixel_info must contain 'pixel_ID'.")
  }

  if (ncol(embedding) < 2) {
    stop("embedding must contain at least 2 columns.")
  }

  out <- data.frame(
    UMAP1 = embedding[[1]],
    UMAP2 = embedding[[2]],
    cluster = cluster,
    run = if ("run" %in% colnames(pixel_info)) pixel_info$run else NA,
    x = if ("x" %in% colnames(pixel_info)) pixel_info$x else NA,
    y = if ("y" %in% colnames(pixel_info)) pixel_info$y else NA,
    z = if ("z" %in% colnames(pixel_info)) pixel_info$z else NA,
    pixel_ID = pixel_info$pixel_ID,
    stringsAsFactors = FALSE
  )

  rownames(out) <- as.character(pixel_info$pixel_ID)

  out
}


#' Attach clustering results back to Cardinal pixelData
#'
#' @param msi_obj A Cardinal MSI object.
#' @param cluster_df A data.frame containing at least:
#'   \code{pixel_ID}, \code{UMAP1}, \code{UMAP2}, \code{kmeans_cluster}.
#'
#' @return The input MSI object with updated pixelData.
#' @export
attach_cluster_to_pixeldata <- function(msi_obj, cluster_df) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }

  pd_df <- as.data.frame(Cardinal::pixelData(msi_obj))

  if (!is.data.frame(cluster_df)) {
    cluster_df <- as.data.frame(cluster_df)
  }

  required_cols <- c("pixel_ID", "UMAP1", "UMAP2", "cluster")
  missing_cols <- setdiff(required_cols, colnames(cluster_df))
  if (length(missing_cols) > 0) {
    stop(
      "cluster_df is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (anyDuplicated(cluster_df$pixel_ID)) {
    stop("cluster_df contains duplicated pixel_ID values.")
  }

  if (!"pixel_ID" %in% colnames(pd_df)) {
    if (is.null(rownames(pd_df))) {
      rownames(pd_df) <- paste0("spectrum=", seq_len(nrow(pd_df)) - 1L)
    }

    pixel_rowname <- rownames(pd_df)

    if ("run" %in% colnames(pd_df)) {
      pixel_ID <- paste(pd_df$run, pixel_rowname, sep = "_")
    } else {
      pixel_ID <- pixel_rowname
    }

    Cardinal::pixelData(msi_obj)$pixel_ID <- pixel_ID
    pd_df <- as.data.frame(Cardinal::pixelData(msi_obj))
  }

  idx <- match(pd_df$pixel_ID, cluster_df$pixel_ID)

  if (anyNA(idx)) {
    stop("Some pixel_ID values in pixelData(msi_obj) are missing from cluster_df.")
  }

  cluster_df2 <- cluster_df[idx, , drop = FALSE]

  if (!all(pd_df$pixel_ID == cluster_df2$pixel_ID)) {
    stop("pixel_ID mismatch after reordering cluster_df.")
  }

  Cardinal::pixelData(msi_obj)$UMAP1 <- cluster_df2$UMAP1
  Cardinal::pixelData(msi_obj)$UMAP2 <- cluster_df2$UMAP2
  Cardinal::pixelData(msi_obj)$cluster <- cluster_df2$cluster

  msi_obj
}

