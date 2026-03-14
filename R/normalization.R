#' Log-transform a matrix
#'
#' Apply log transformation to stabilize variance and compress the dynamic
#' range of intensity values. The transformation used is:
#'
#' \deqn{log10(x + 1)}
#'
#' This step is commonly applied after RMS or TIC normalization in
#' mass spectrometry imaging workflows.
#'
#' @param x A numeric matrix with pixels/samples as rows and features as columns.
#'
#' @return A numeric matrix of the same dimension as \code{x}.
#'
#' @export
log_transform <- function(x) {
  check_spectra_matrix(x)
  log10(x + 1)
}


#' Z-score scale a matrix by columns
#'
#' Perform column-wise z-score scaling:
#'
#' \deqn{(x - mean(x)) / sd(x)}
#'
#' Each feature is centered to mean 0 and scaled to unit variance.
#'
#' @param x A numeric matrix with pixels/samples as rows and features as columns.
#'
#' @return A numeric matrix where each column has mean 0 and standard deviation 1.
#'
#' @export
zscore_scale <- function(x) {

  check_spectra_matrix(x)

  means <- colMeans(x, na.rm = TRUE)
  sds <- apply(x, 2, stats::sd, na.rm = TRUE)

  if (any(sds == 0)) {
    stop(
      "Matrix contains constant columns. ",
      "Please remove them using remove_constant_features()."
    )
  }

  x_scaled <- sweep(sweep(x, 2, means, "-"), 2, sds, "/")
  x_scaled <- as.matrix(x_scaled)

  rownames(x_scaled) <- rownames(x)
  colnames(x_scaled) <- colnames(x)

  x_scaled
}


#' Min-max scale a matrix by columns
#'
#' Perform column-wise min-max scaling:
#'
#' \deqn{(x - min(x)) / (max(x) - min(x))}
#'
#' Each feature is rescaled independently to the range \code{[0,1]}.
#'
#' @param x A numeric matrix with pixels/samples as rows and features as columns.
#'
#' @return A numeric matrix of the same dimension as \code{x}.
#'
#' @export
minmax_scale <- function(x) {

  check_spectra_matrix(x)

  min_vals <- apply(x, 2, min, na.rm = TRUE)
  max_vals <- apply(x, 2, max, na.rm = TRUE)
  ranges <- max_vals - min_vals

  if (any(ranges == 0)) {
    stop(
      "Matrix contains constant columns. ",
      "Please remove them using remove_constant_features()."
    )
  }

  x_scaled <- sweep(sweep(x, 2, min_vals, "-"), 2, ranges, "/")
  x_scaled <- as.matrix(x_scaled)

  rownames(x_scaled) <- rownames(x)
  colnames(x_scaled) <- colnames(x)

  x_scaled
}


#' Apply feature scaling or transformation
#'
#' Apply optional feature scaling or transformation to a matrix after
#' spectral normalization (e.g., RMS or TIC normalization).
#'
#' Supported methods include:
#'
#' \itemize{
#'   \item \code{"log"}: log10(x + 1) transformation
#'   \item \code{"zscore"}: column-wise z-score scaling
#'   \item \code{"minmax"}: column-wise min-max scaling
#'   \item \code{"none"}: no transformation
#' }
#'
#' @param x A numeric matrix with pixels/samples as rows and features as columns.
#' @param method Character string specifying the scaling method.
#'
#' @return A numeric matrix of the same dimension as \code{x}.
#'
#' @examples
#' mat <- matrix(runif(100), nrow = 10)
#'
#' apply_feature_scaling(mat, "log")
#' apply_feature_scaling(mat, "zscore")
#' apply_feature_scaling(mat, "minmax")
#'
#' @export
apply_feature_scaling <- function(x, method = c("log", "zscore", "minmax", "none")) {

  check_spectra_matrix(x)

  method <- match.arg(method)

  if (method == "log") {
    return(log_transform(x))
  }

  if (method == "zscore") {
    return(zscore_scale(x))
  }

  if (method == "minmax") {
    return(minmax_scale(x))
  }

  if (method == "none") {
    return(x)
  }
}

