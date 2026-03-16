#' Cherry tomato MSI example data
#'
#' A toy `MSImagingExperiment` dataset included in \pkg{spatMetaCluster} for
#' demonstrating spatial metabolomics workflows.
#'
#' This dataset is manually constructed for demonstration purposes and is not
#' derived from a real mass spectrometry imaging experiment. Its spatial pattern
#' and RGB pixel information were generated based on a real cherry tomato image,
#' while the m/z features were manually defined for illustrative use only.
#'
#' Accordingly, the m/z values in this dataset do not represent real metabolites
#' and should not be interpreted as having chemical or biological meaning. The
#' dataset is intended exclusively for examples, tutorials, testing, and package
#' demonstrations.
#'
#' @format An object of class `MSImagingExperiment` with 32 features and 43,253
#' spectra.
#'
#' \describe{
#'   \item{spectraData}{One intensity assay.}
#'   \item{featureData}{Feature metadata including `mz`, `count`, and `freq`.}
#'   \item{pixelData}{Pixel metadata including spatial coordinates (`x`, `y`),
#'   run information, cluster information(`Fruit`, `Calyx`), and `pixel_ID`.}
#'   \item{mass range}{From 140.3475 to 259.6539.}
#'   \item{centroided}{`TRUE`.}
#' }
#'
#' @details
#' This lightweight dataset is designed to provide a reproducible example for
#' data handling, visualization, and clustering workflows implemented in
#' \pkg{spatMetaCluster}.
#'
#' @source Constructed by the package author from a self-acquired cherry tomato
#' photograph, with manually defined example m/z features for demonstration only.
#'
#' @usage data(cherry_tomato_msi)
#' @docType data
#' @keywords datasets
#' @name cherry_tomato_msi
NULL
