#' scOverlay: layered overlays for single-cell embeddings
#'
#' scOverlay creates layered overlay plots for reduced-dimension embeddings
#' stored in \linkS4class{SingleCellExperiment} objects. It combines background
#' and foreground layers to show gene expression, assay-derived features and
#' \code{colData} metadata with independent colour scales.
#'
#' The package includes helpers for single plots, grouped plots and nested
#' grouped plots; bivariate colour blends for two numeric foreground variables;
#' expression-combination plots for small sets of genes or numeric metadata; and
#' palette listing and preview utilities. Background and foreground point layers
#' can be rasterized independently for dense single-cell plots while keeping
#' labels, axes and legends as vector graphics.
#'
#' @name scOverlay
#' @keywords internal
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment rownames
#' @importFrom ggplot2 ggplot scale_colour_gradientn scale_colour_manual
#' @importFrom ggplot2 theme_classic theme_void labs annotate aes
#' @importFrom rlang .data
#' @importFrom stats na.omit

NULL
