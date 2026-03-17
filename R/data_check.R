#' Extract spectra matrix and pixel metadata from a Cardinal MSI object
#'
#' @param msi_obj A Cardinal MSI object imported by \code{Cardinal::readImzML()}.
#'
#' @return A list containing:
#' \describe{
#'   \item{spectra}{Numeric matrix with pixels as rows and m/z features as columns.}
#'   \item{pixel_info}{A data.frame of pixel metadata.}
#'   \item{mz}{Numeric vector of m/z values.}
#' }
#'
#' @export
extract_spectra_matrix <- function(msi_obj) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }

  intensity_obj <- Cardinal::intensity(msi_obj)
  pixel_info <- as.data.frame(Cardinal::pixelData(msi_obj))
  mz_vals <- Cardinal::mz(msi_obj)

  nr <- nrow(intensity_obj)
  nc <- ncol(intensity_obj)

  intensity_mat <- tryCatch(
    {
      tmp <- intensity_obj[seq_len(nr), seq_len(nc)]
      as.matrix(tmp)
    },
    error = function(e) {
      stop(
        "Failed to extract intensity matrix from Cardinal object. ",
        "Original error: ", conditionMessage(e)
      )
    }
  )

  if (!is.matrix(intensity_mat)) {
    stop("Intensity data could not be converted to a standard matrix.")
  }

  spectra_mat <- t(intensity_mat)
  storage.mode(spectra_mat) <- "numeric"

  if (ncol(spectra_mat) != length(mz_vals)) {
    stop(
      "Number of columns in spectra matrix does not match length of mz values: ",
      ncol(spectra_mat), " vs ", length(mz_vals), "."
    )
  }

  colnames(spectra_mat) <- as.character(mz_vals)

  if (nrow(spectra_mat) != nrow(pixel_info)) {
    stop(
      "Number of rows in spectra matrix does not match rows in pixelData: ",
      nrow(spectra_mat), " vs ", nrow(pixel_info), "."
    )
  }

  if (!"pixel_ID" %in% colnames(pixel_info)) {
    if (is.null(rownames(pixel_info))) {
      rownames(pixel_info) <- paste0("spectrum=", seq_len(nrow(pixel_info)) - 1L)
    }

    pixel_info$pixel_rowname <- rownames(pixel_info)

    if ("run" %in% colnames(pixel_info)) {
      pixel_info$pixel_ID <- paste(pixel_info$run, seq_len(nrow(pixel_info)), sep = "_")
    } else {
      pixel_info$pixel_ID <- paste0("pixel_", seq_len(nrow(pixel_info)))
    }
  }

  rownames(spectra_mat) <- pixel_info$pixel_ID

  list(
    spectra = spectra_mat,
    pixel_info = pixel_info,
    mz = mz_vals
  )
}


#' Check whether a spectra matrix is valid
#'
#' This function performs basic checks on a spectra matrix before downstream
#' preprocessing or clustering.
#'
#' @param x A numeric matrix with samples/pixels as rows and features as columns.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#'
#' @examples
#' mat <- matrix(runif(20), nrow = 5)
#' check_spectra_matrix(mat)
#'
#' @export
check_spectra_matrix <- function(x) {
  if (!is.matrix(x)) {
    stop("Input 'x' must be a matrix.")
  }

  if (!is.numeric(x)) {
    stop("Input matrix must be numeric.")
  }

  if (nrow(x) < 2) {
    stop("Input matrix must have at least 2 rows.")
  }

  if (ncol(x) < 2) {
    stop("Input matrix must have at least 2 columns.")
  }

  if (any(is.na(x))) {
    warning("Input matrix contains NA values.")
  }

  invisible(TRUE)
}

#' Remove constant features from a spectra matrix
#'
#' This function removes columns with zero standard deviation, which can cause
#' issues in downstream normalization or dimensionality reduction.
#'
#' @param x A numeric matrix with samples/pixels as rows and features as columns.
#'
#' @return A numeric matrix with constant columns removed.
#'
#' @examples
#' mat <- cbind(
#'   a = c(1, 2, 3),
#'   b = c(5, 5, 5),
#'   c = c(2, 4, 6)
#' )
#' remove_constant_features(mat)
#'
#' @export
remove_constant_features <- function(x) {
  check_spectra_matrix(x)

  keep <- apply(x, 2, stats::sd, na.rm = TRUE) > 0
  x_filtered <- x[, keep, drop = FALSE]

  if (ncol(x_filtered) == 0) {
    stop("All columns were removed because all features are constant.")
  }

  x_filtered
}
