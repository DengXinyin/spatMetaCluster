#' Extract spectra matrix and pixel metadata from a Cardinal MSI object
#'
#' @param msi_obj A Cardinal MSI object imported by `Cardinal::readImzML()`.
#'
#' @return A list containing:
#' \describe{
#'   \item{spectra}{Numeric matrix with pixels as rows and m/z features as columns.}
#'   \item{pixel_info}{A data.frame of pixel metadata.}
#'   \item{mz}{Numeric vector of m/z values.}
#'   \item{msi_obj}{The MSI object with ensured `pixel_ID` in `pixelData`.}
#' }
#'
#' @export
extract_spectra_matrix <- function(msi_obj) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }

  msi_obj <- ensure_pixel_id(msi_obj)

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
    stop("Internal error: 'pixel_ID' should exist after ensure_pixel_id().")
  }

  rownames(spectra_mat) <- pixel_info$pixel_ID

  list(
    spectra = spectra_mat,
    pixel_info = pixel_info,
    mz = mz_vals,
    msi_obj = msi_obj
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



#' Ensure a Cardinal MSI object contains a pixel_ID column
#'
#' If `pixelData(msi_obj)` does not contain a `pixel_ID` column, this function
#' creates one and writes it back into the MSI object.
#'
#' @param msi_obj A Cardinal MSI object.
#'
#' @return The same MSI object with `pixel_ID` present in `pixelData`.
#' @export
ensure_pixel_id <- function(msi_obj) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }

  pd <- Cardinal::pixelData(msi_obj)
  pd_df <- as.data.frame(pd)

  if (!"pixel_ID" %in% colnames(pd_df)) {
    if ("run" %in% colnames(pd_df)) {
      pd$pixel_ID <- paste(pd_df$run, seq_len(nrow(pd_df)), sep = "_")
    } else {
      pd$pixel_ID <- paste0("pixel_", seq_len(nrow(pd_df)))
    }

    Cardinal::pixelData(msi_obj) <- pd
  }

  msi_obj
}
