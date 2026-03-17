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
#'   Default is \code{median}.
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
#' @details
#' Spatially enriched metabolites are defined as features satisfying both
#' high specificity and high abundance:
#' \itemize{
#'   \item \code{is_high_specific}: \code{|statistic| > t_statistic_threshold},
#'   \code{|cor| > cor_threshold}, and \code{M2 > M2_threshold}
#'   \item \code{is_high_abundance}: \code{M1 > M1_threshold}
#' }
#'
#' The final filtered table is:
#' \preformatted{
#' SEMs_df <- results_df[
#'   results_df$is_high_specific | results_df$is_high_abundance,
#'   , drop = FALSE
#' ]
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' data(cherry_tomato_msi)
#' res <- SEMs_screen(
#'   msi_obj = cherry_tomato_msi,
#'   group_col = "ROI"
#' )
#' head(res$SEMs_df)
#' }
#'
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
    coloc_threshold = median
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

  group_vec <- pixel_df[[group_col]]

  if (all(is.na(group_vec))) {
    stop("Grouping column '", group_col, "' contains only NA values.")
  }

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

  coloc_names <- names(coloc_res)

  colocalized_df <- dplyr::bind_rows(lapply(coloc_names, function(nm) {
    df <- as.data.frame(coloc_res[[nm]])
    df$tissue <- nm
    df
  }))

  merged_df <- dplyr::inner_join(
    ssc_df,
    colocalized_df,
    by = c("i", "mz")
  )

  if (!all(c("class", "tissue") %in% colnames(merged_df))) {
    stop("Merged table must contain both 'class' and 'tissue' columns.")
  }

  matched_df <- merged_df[merged_df$class == merged_df$tissue, , drop = FALSE]

  results_df <- matched_df %>%
    dplyr::mutate(
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
    , drop = FALSE
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

