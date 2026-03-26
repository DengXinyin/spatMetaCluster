#' Run UMAP using a Python backend via reticulate
#'
#' This function calls Python \pkg{umap-learn} through \pkg{reticulate} and
#' returns the embedding as an R data.frame.
#'
#' A working Python environment with \pkg{umap-learn} and \pkg{numpy} installed
#' is required. This function was tested with Python 3.9, \pkg{numpy} 1.24.4,
#' and \pkg{umap-learn} 0.5.7 in a conda environment, although other compatible
#' versions may also work.
#'
#' @param x A numeric matrix with samples/pixels as rows and features as columns.
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#' @param metric Distance metric used by UMAP. Must be one of
#' \code{"cosine"} or \code{"euclidean"}. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors used by UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of embedding dimensions. Default is \code{2L}.
#' Any positive integer is allowed, but \code{2} or \code{3} is recommended
#' for most visualization and exploratory analysis tasks.
#' @param random_state Optional integer; random seed for UMAP. Default is
#' \code{2025L}. Use \code{NULL} to allow non-deterministic execution. In the
#' Python \pkg{umap-learn} backend, setting a fixed \code{random_state} may
#' disable or limit parallel execution.
#' @param n_jobs Integer; number of parallel jobs used by UMAP. Default is
#' \code{1L}. Use \code{-1L} to request all available CPU cores if supported
#' by the Python backend. If \code{n_jobs != 1L} and \code{random_state} is not
#' \code{NULL}, this function will issue a warning and set the Python UMAP
#' \code{random_state} to \code{NULL} to allow parallel execution.
#' @param verbose Logical; whether to print UMAP progress. Default is \code{TRUE}.
#'
#' @return A data.frame containing UMAP coordinates. Column names are
#' \code{UMAP1}, \code{UMAP2}, ... depending on \code{n_components}.
#'
#' @details
#' In Python \pkg{umap-learn}, fixed random seeds and parallel execution may be
#' incompatible. In particular, when \code{random_state} is fixed, the backend
#' may override \code{n_jobs} to \code{1}. To enable actual parallel execution,
#' use \code{random_state = NULL}.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(runif(100), nrow = 10)
#'
#' emb_single <- run_umap_py(
#'   mat,
#'   python_path = "/path/to/python",
#'   n_components = 2,
#'   random_state = 2025L,
#'   n_jobs = 1L
#' )
#' head(emb_single)
#'
#' emb_parallel <- run_umap_py(
#'   mat,
#'   python_path = "/path/to/python",
#'   metric = "euclidean",
#'   n_components = 2,
#'   random_state = NULL,
#'   n_jobs = 4L
#' )
#' head(emb_parallel)
#' }
#'
#' @export
run_umap_py <- function(
    x,
    python_path = NULL,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    random_state = 2025L,
    n_jobs = 1L,
    verbose = TRUE
) {
  check_spectra_matrix(x)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  metric <- match.arg(metric, choices = c("cosine", "euclidean"))

  if (!is.numeric(n_neighbors) || length(n_neighbors) != 1L || is.na(n_neighbors)) {
    stop("'n_neighbors' must be a single positive integer.")
  }
  n_neighbors <- as.integer(n_neighbors)
  if (n_neighbors < 2L) {
    stop("'n_neighbors' must be >= 2.")
  }

  if (!is.numeric(min_dist) || length(min_dist) != 1L || is.na(min_dist)) {
    stop("'min_dist' must be a single non-negative numeric value.")
  }
  if (min_dist < 0) {
    stop("'min_dist' must be >= 0.")
  }

  if (!is.numeric(n_components) || length(n_components) != 1L || is.na(n_components)) {
    stop("'n_components' must be a single positive integer.")
  }
  n_components <- as.integer(n_components)
  if (n_components < 1L) {
    stop("'n_components' must be a positive integer.")
  }

  if (!is.null(random_state)) {
    if (!is.numeric(random_state) || length(random_state) != 1L || is.na(random_state)) {
      stop("'random_state' must be NULL or a single integer.")
    }
    random_state <- as.integer(random_state)
  }

  if (!is.numeric(n_jobs) || length(n_jobs) != 1L || is.na(n_jobs)) {
    stop("'n_jobs' must be a single integer.")
  }
  n_jobs <- as.integer(n_jobs)
  if (n_jobs == 0L) {
    stop("'n_jobs' must not be 0. Use 1L for single-threading or -1L for all available cores if supported.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  if (n_jobs != 1L && !is.null(random_state)) {
    warning(
      "In Python 'umap-learn', a fixed 'random_state' may force 'n_jobs' to 1. ",
      "To allow parallel execution, 'random_state' will be set to NULL for the Python UMAP call."
    )
    random_state_py <- NULL
  } else {
    random_state_py <- random_state
  }

  umap <- reticulate::import("umap", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  if (!is.null(random_state_py)) {
    np$random$seed(random_state_py)
  }

  umap_result <- umap$UMAP(
    metric = metric,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    random_state = random_state_py,
    n_jobs = n_jobs,
    verbose = verbose
  )$fit_transform(x)

  umap_matrix <- reticulate::py_to_r(
    np$array(umap_result, dtype = np$float64)
  )

  umap_df <- as.data.frame(umap_matrix)
  colnames(umap_df) <- paste0("UMAP", seq_len(n_components))
  rownames(umap_df) <- rownames(x)

  umap_df
}


#' Check Python UMAP availability
#'
#' This helper function checks whether the required Python packages
#' \pkg{umap-learn} and \pkg{numpy} are available in the current or specified
#' Python environment.
#'
#' The Python backend was tested with Python 3.9, \pkg{numpy} 1.24.4, and
#' \pkg{umap-learn} 0.5.7 in a conda environment, although other compatible
#' versions may also work.
#'
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#'
#' @return Invisibly returns \code{TRUE} if Python UMAP dependencies are available.
#'
#' @examples
#' \dontrun{
#' check_umap_python_env("/path/to/python")
#' }
#'
#' @export
check_umap_python_env <- function(python_path = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  py_config <- reticulate::py_config()

  if (!reticulate::py_module_available("umap")) {
    stop("Python module 'umap' is not available in the selected environment.")
  }

  if (!reticulate::py_module_available("numpy")) {
    stop("Python module 'numpy' is not available in the selected environment.")
  }

  message("Python environment OK: ", py_config$python)
  invisible(TRUE)
}
