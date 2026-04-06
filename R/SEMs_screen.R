#' Screen spatially enriched metabolites using SSC and colocalization
#'
#' This function identifies spatially enriched metabolites (SEMs) by combining
#' spatial shrunken centroids (SSC) results and spatial colocalization analysis
#' from a Cardinal MSI object. Only the final SEMs table is saved when
#' \code{save_table = TRUE}.
#'
#' @param msi_obj A Cardinal MSI object.
#' @param group_col A character string specifying the column name in
#'   \code{pixelData(msi_obj)} used as the grouping variable, such as \code{"ROI"}.
#' @param t_statistic_threshold Numeric; threshold for absolute SSC t statistic.
#'   Default is \code{20}.
#' @param cor_threshold Numeric; threshold for absolute colocalization
#'   correlation. Default is \code{0.1}.
#' @param M1_threshold Numeric; threshold for M1, used to define abundance.
#'   Default is \code{0.2}.
#' @param M2_threshold Numeric; threshold for M2, used to define specificity.
#'   Default is \code{0.5}.
#' @param save_table Logical; whether to save the final \code{SEMs_df}.
#'   Default is \code{FALSE}.
#' @param file_name Character string specifying the output file name without
#'   extension. If \code{NULL}, the file name will be automatically generated as
#'   \code{"Spatially_Enriched_Metabolites_of_<msi_object_name>"}.
#' @param file_format Character string specifying output file format. Must be
#'   one of \code{"csv"} or \code{"xlsx"}. Default is \code{"csv"}.
#' @param ssc_seed Integer; random seed for SSC analysis. Default is \code{2025L}.
#' @param coloc_threshold Threshold function passed to \code{colocalized()}.
#'   Default is \code{stats::median}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{ssc_df}{A data.frame of SSC top features.}
#'   \item{colocalized_df}{A data.frame of colocalization results.}
#'   \item{merged_df}{Merged SSC and colocalization table.}
#'   \item{matched_df}{Rows where SSC class matches colocalization tissue label.}
#'   \item{results_df}{Annotated screening table before SEM filtering.}
#'   \item{SEMs_df}{Final filtered table of spatially enriched metabolites.}
#' }
#'
#' @importFrom rlang .data
#' @export
SEMs_screen <- function(
    msi_obj,
    group_col,
    t_statistic_threshold = 20,
    cor_threshold = 0.1,
    M1_threshold = 0.2,
    M2_threshold = 0.5,
    save_table = FALSE,
    file_name = NULL,
    file_format = "csv",
    ssc_seed = 2025L,
    coloc_threshold = stats::median
) {
  if (!is.character(group_col) || length(group_col) != 1L) {
    stop("'group_col' must be a single character string.")
  }

  if (!file_format %in% c("csv", "xlsx")) {
    stop("'file_format' must be either 'csv' or 'xlsx'.")
  }

  msi_name <- deparse(substitute(msi_obj))
  if (!is.character(msi_name) || length(msi_name) != 1L || nchar(msi_name) == 0L) {
    msi_name <- "MSI_data"
  }

  if (is.null(file_name)) {
    file_name <- paste0("Spatially_Enriched_Metabolites_of_", msi_name)
  }
  file_name <- gsub("[^A-Za-z0-9_\\-]", "_", file_name)

  pixel_df <- as.data.frame(Cardinal::pixelData(msi_obj))

  if (!group_col %in% colnames(pixel_df)) {
    stop("Column '", group_col, "' was not found in pixelData(msi_obj).")
  }

  group_raw <- pixel_df[[group_col]]

  if (all(is.na(group_raw))) {
    stop("Grouping column '", group_col, "' contains only NA values.")
  }

  non_na_groups <- unique(group_raw[!is.na(group_raw)])
  if (length(non_na_groups) < 2) {
    stop("Grouping column '", group_col, "' must contain at least 2 non-NA groups.")
  }

  group_vec <- as.factor(group_raw)

  set.seed(ssc_seed)
  ssc_res <- Cardinal::spatialShrunkenCentroids(
    msi_obj,
    y = group_vec,
    weights = "adaptive",
    verbose = Cardinal::getCardinalVerbose()
  )

  ssc_df <- as.data.frame(Cardinal::topFeatures(
    ssc_res,
    sort.by = "statistic"
  ))

  coloc_res <- Cardinal::colocalized(
    msi_obj,
    ref = group_vec,
    threshold = coloc_threshold,
    n = Inf,
    sort.by = "cor",
    verbose = Cardinal::getCardinalVerbose(),
    chunkopts = list(),
    BPPARAM = Cardinal::getCardinalBPPARAM()
  )

  colocalized_df <- NULL

  if (is.list(coloc_res) &&
      !inherits(coloc_res, "DataFrame") &&
      !inherits(coloc_res, "data.frame")) {

    coloc_names <- names(coloc_res)

    if (is.null(coloc_names) || any(coloc_names == "")) {
      stop(
        "Cardinal::colocalized() returned a list, but list elements are not named. ",
        "Unable to determine tissue/group labels."
      )
    }

    colocalized_df <- dplyr::bind_rows(lapply(seq_along(coloc_res), function(i) {
      df <- as.data.frame(coloc_res[[i]])
      df$tissue <- coloc_names[i]
      df
    }))
  } else {
    colocalized_df <- as.data.frame(coloc_res)

    if (!"tissue" %in% colnames(colocalized_df)) {
      if ("class" %in% colnames(colocalized_df)) {
        colocalized_df$tissue <- as.character(colocalized_df$class)
      } else if ("group" %in% colnames(colocalized_df)) {
        colocalized_df$tissue <- as.character(colocalized_df$group)
      } else if ("ref" %in% colnames(colocalized_df)) {
        colocalized_df$tissue <- as.character(colocalized_df$ref)
      } else {
        stop(
          "Could not determine group labels from the output of Cardinal::colocalized(). ",
          "Expected either:\n",
          "1) a named list of group-wise results, or\n",
          "2) a data.frame/DataFrame containing a 'tissue', 'class', 'group', or 'ref' column.\n",
          "Actual columns were: ",
          paste(colnames(colocalized_df), collapse = ", ")
        )
      }
    }
  }

  required_ssc_cols <- c("i", "mz", "class", "statistic")
  missing_ssc_cols <- setdiff(required_ssc_cols, colnames(ssc_df))
  if (length(missing_ssc_cols) > 0) {
    stop(
      "SSC result is missing required columns: ",
      paste(missing_ssc_cols, collapse = ", ")
    )
  }

  required_coloc_cols <- c("i", "mz", "cor", "M1", "M2", "tissue")
  missing_coloc_cols <- setdiff(required_coloc_cols, colnames(colocalized_df))
  if (length(missing_coloc_cols) > 0) {
    stop(
      "Colocalization result is missing required columns: ",
      paste(missing_coloc_cols, collapse = ", ")
    )
  }

  merged_df <- dplyr::inner_join(
    ssc_df,
    colocalized_df,
    by = c("i", "mz")
  )

  if (!all(c("class", "tissue") %in% colnames(merged_df))) {
    stop("Merged table must contain both 'class' and 'tissue' columns.")
  }

  matched_df <- merged_df[
    as.character(merged_df$class) == as.character(merged_df$tissue),
    ,
    drop = FALSE
  ]

  results_df <- dplyr::mutate(
    matched_df,
    is_high_specific =
      (abs(.data$statistic) > t_statistic_threshold) &
      (abs(.data$cor) > cor_threshold) &
      (.data$M2 > M2_threshold),

    is_high_abundance =
      (.data$M1 > M1_threshold),

    specificity_tag = dplyr::case_when(
      .data$statistic > t_statistic_threshold &
        .data$cor > cor_threshold &
        .data$M2 > M2_threshold ~ "Specifically enriched",

      .data$statistic < -t_statistic_threshold &
        .data$cor < -cor_threshold ~ "Specifically depleted",

      abs(.data$statistic) < t_statistic_threshold / 2 &
        .data$M1 < M1_threshold / 2 ~ "No clear spatial pattern",

      TRUE ~ "Complex distribution"
    ),

    abundance_tag = dplyr::if_else(
      .data$M1 > M1_threshold,
      "High abundance",
      "Low abundance"
    ),

    composite_tag = paste(.data$specificity_tag, .data$abundance_tag, sep = " | ")
  )

  SEMs_df <- results_df[
    results_df$is_high_specific & results_df$is_high_abundance,
    ,
    drop = FALSE
  ]

  if (save_table) {
    out_file <- paste0(file_name, ".", file_format)

    if (file_format == "csv") {
      utils::write.csv(SEMs_df, file = out_file, row.names = FALSE)
    }

    if (file_format == "xlsx") {
      if (!requireNamespace("writexl", quietly = TRUE)) {
        stop(
          "Package 'writexl' is required to save xlsx files. ",
          "Please install it first."
        )
      }
      writexl::write_xlsx(SEMs_df, path = out_file)
    }
  }

  list(
    ssc_df = ssc_df,
    colocalized_df = colocalized_df,
    merged_df = merged_df,
    matched_df = matched_df,
    results_df = results_df,
    SEMs_df = SEMs_df
  )
}
