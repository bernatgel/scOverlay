#' Small example SingleCellExperiment for scOverlay
#'
#' A small \linkS4class[SingleCellExperiment]{SingleCellExperiment} object used in scOverlay examples
#' and tests. It contains a small synthetic set of cells, genes, cell
#' annotations and reduced dimensions to demonstrate the plotting functions
#' quickly.
#' It includes \code{logcounts}, \code{TSNE} and \code{UMAP} reduced dimensions,
#' and cell metadata such as \code{sample}, \code{sample_type} and
#' \code{cluster}.
#'
#' @format A \linkS4class[SingleCellExperiment]{SingleCellExperiment} object with gene expression
#' values, cell metadata and reduced dimensions.
#'
#' @examples
#' data("sce_overlay_example")
#' sce_overlay_example
#'
#' @name sce_overlay_example
#' @docType data
#' @keywords datasets
NULL

#' Example MPNST SingleCellExperiment for the scOverlay vignette
#'
#' A small \linkS4class[SingleCellExperiment]{SingleCellExperiment} object derived from the single-cell
#' MPNST progression dataset deposited at the European Genome-phenome Archive
#' under accession EGAS50000001747. The object has been reduced in size and is
#' included only to demonstrate scOverlay functionality in the package vignette.
#'
#' @format A \linkS4class[SingleCellExperiment]{SingleCellExperiment} object with selected genes, cell
#' metadata and reduced dimensions.
#'
#' @source European Genome-phenome Archive study EGAS50000001747,
#' \url{https://ega-archive.org/studies/EGAS50000001747}.
#'
#' @examples
#' data("sce_mpnst_example")
#' sce_mpnst_example
#'
#' @name sce_mpnst_example
#' @docType data
#' @keywords datasets
NULL
