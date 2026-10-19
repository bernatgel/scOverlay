#' Plot overlays for several genes
#'
#' \code{plotGeneOverlay} creates one \code{\link{plotOverlay}} plot for each
#' gene in \code{genes} and returns them as a named list.
#'
#' This is a convenience function to quickly inspect the expression pattern of
#' several genes on the same embedding. All additional plotting arguments are
#' passed to \code{\link{plotOverlay}}.
#'
#' @details
#' Genes not present in \code{rownames(sce)} are skipped with a warning. If none
#' of the requested genes are found, the function stops with an error.
#'
#' @param sce A \linkS4class[SingleCellExperiment]{SingleCellExperiment} object containing the reduced
#' dimensions and gene expression values to plot.
#'
#' @param genes Character vector with the genes to plot. Genes are looked up in
#' \code{rownames(sce)}. This argument has no default and must be provided by the
#' user.
#'
#' @param ... Additional arguments passed to \code{\link{plotOverlay}}, such as
#' \code{reduced_dim}, \code{background}, \code{bg_type}, \code{bg_assay},
#' \code{fg_assay}, \code{fg_order}, palettes or rasterization options.
#'
#' @return A named list of \pkg{ggplot2} objects, one for each gene found in
#' \code{rownames(sce)}.
#'
#' @seealso \code{\link{plotOverlay}}, \code{\link{getPalette}} and
#' \code{\link{listPalettes}}.
#'
#' @examples
#' data("sce_overlay_example")
#'
#' plots <- plotGeneOverlay(
#'     sce = sce_overlay_example,
#'     genes = c("SOX10", "S100B"),
#'     reduced_dim = "TSNE",
#'     background = "cluster",
#'     fg_order = "ascending"
#' )
#'
#' plots$SOX10
#'
#' @export



plotGeneOverlay <- function(sce, genes, ...) {
	dots <- list(...)
	dot_names <- names(dots)
	if (is.null(dot_names)) {
		dot_names <- character(0L)
	}
	if ("foreground" %in% dot_names) {
		stop(
			"`foreground` cannot be supplied through `...` in `plotGeneOverlay()`; use `genes` instead.",
			call. = FALSE
		)
	}

	.validate_sce(sce)
	if (!is.character(genes) || length(genes) == 0L || anyNA(genes)) {
		stop("`genes` must be a non-empty character vector without NA values.", call. = FALSE)
	}

	available_genes <- genes[genes %in% rownames(sce)]
	if (length(available_genes) == 0L) {
		stop("None of the requested `genes` were found in rownames(sce).", call. = FALSE)
	}

	missing_genes <- setdiff(genes, rownames(sce))
	if (length(missing_genes) > 0L) {
		warning(
			"`genes` contains features not present in rownames(sce) and they will be skipped: ",
			paste(missing_genes, collapse = ", "),
			".",
			call. = FALSE
		)
	}

	genes <- available_genes
	
	plots <- lapply(genes, function(g) {
		return(plotOverlay(sce = sce, foreground = g, ...))
	})
	names(plots) <- genes
	return(plots)
}
