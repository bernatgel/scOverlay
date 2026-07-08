#' Plot overlays split by a cell metadata group
#'
#' \code{plotOverlayPerGroup} creates one \code{\link{plotOverlay}} plot for
#' each group defined by a column in \code{colData(sce)} and combines the plots
#' into a single \pkg{patchwork} object.
#'
#' This function is useful to compare the same foreground variable across
#' samples, sample types, clusters or any other cell-level annotation. Each
#' panel keeps the full embedding as background, while the foreground layer is
#' restricted to the cells belonging to the corresponding group.
#'
#' @details
#' The grouping variable is specified with \code{group_col}. For each value in
#' this column, \code{plotOverlayPerGroup} calls \code{\link{plotOverlay}} using
#' the same background, foreground, palette and plotting options. The foreground
#' subset for each panel is defined by the group membership and, optionally, by
#' the additional filter supplied with \code{fg_subset_cells}.
#'
#' As in \code{\link{plotOverlay}}, subsetting only affects the foreground
#' layer. Unless \code{background = "none"}, the background layer still shows all
#' cells in every panel, providing a common reference embedding across groups.
#' This makes it possible to compare where cells from different samples or
#' sample types are located in the same reduced dimension space.
#'
#' With \code{fg_limits = "shared"}, foreground ranges are computed after the
#' global \code{fg_subset_cells} filter and before splitting into groups. In
#' bivariate mode, separate shared ranges are computed for the two numeric
#' foreground variables so colours are comparable across panels. With
#' \code{bg_limits = "shared"}, one background range is computed across all
#' cells. Categorical foreground and background palettes are based on global
#' categorical values, so colours remain consistent across panels. Named
#' categorical palettes are recommended when strict category-colour control is
#' needed.
#'
#' If \code{shared_legend = TRUE}, compatible legends are collected with
#' \pkg{patchwork}. The legend can be repositioned afterwards with patchwork
#' syntax, for example \code{p & ggplot2::theme(legend.position = "bottom")}.
#'
#' If \code{group_order} is provided, panels are drawn following that order. If
#' it is \code{NULL}, factor levels are used when \code{group_col} is a factor;
#' otherwise, groups are shown in their order of appearance in the object.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object containing the reduced
#' dimensions and the data to plot. The object must contain the reduced
#' dimension specified in \code{reduced_dim}, the grouping column specified in
#' \code{group_col}, and the data required by \code{foreground}.
#'
#' @param group_col A single character value with the name of the
#' \code{colData(sce)} column used to split the plot into panels. Each value in
#' this column defines one panel. This argument has no default and must be
#' provided by the user.
#'
#' @param foreground Foreground layer specification. It can be a feature name in
#' \code{rownames(sce)}, a \code{colData(sce)} column, \code{"solid"},
#' \code{"none"}, or a character vector of feature or column names for
#' \code{fg_type = "bivariate"} or \code{fg_type = "combination"}. This
#' argument has no default and must be provided by the user.
#'
#' @param background A single character value specifying the background layer
#' passed to \code{\link{plotOverlay}}. It can be \code{"solid"} to draw all
#' background cells with a fixed colour, \code{"none"} to draw no visible
#' background points, a feature in \code{rownames(sce)}, or a column in
#' \code{colData(sce)}. Defaults to \code{"solid"}.
#'
#' @param reduced_dim A single character value with the name of the reduced
#' dimension to plot. It must be one of \code{reducedDimNames(sce)}. Only the
#' first two dimensions are used. Defaults to \code{"TSNE"}.
#'
#' @param group_order Optional character vector specifying the order in which
#' groups should be plotted. If \code{NULL}, factor levels are used for factor
#' grouping variables; otherwise, the order of appearance in \code{group_col} is
#' used. Defaults to \code{NULL}.
#'
#' @param fg_subset_cells An optional additional cell subset applied to the
#' foreground layer before splitting by group. It can be \code{NULL}, a logical
#' vector of length \code{ncol(sce)}, a character vector with cell names, a
#' function receiving \code{sce} and returning a logical vector, or an expression
#' evaluated in \code{colData(sce)}. \code{NULL} selects all cells. \code{NA}
#' values are treated as \code{FALSE}. Defaults to \code{NULL}.
#' Groups with no matching foreground cells still produce panels with empty
#' foreground layers.
#'
#' @param bg_type A character value indicating how the background variable should
#' be represented when \code{background} is not \code{"solid"} or
#' \code{"none"}. Accepted values are \code{"auto"}, \code{"continuous"} and
#' \code{"categorical"}. See \code{\link{plotOverlay}} for details. Defaults to
#' \code{"auto"}.
#'
#' @param bg_assay A single character value with the name of the assay to use
#' when \code{background} is a feature in \code{rownames(sce)}. It is ignored
#' when \code{background = "solid"} or \code{background = "none"}. Defaults to
#' \code{"logcounts"}.
#'
#' @param fg_assay A single character value with the name of the assay to use
#' when \code{foreground} is a feature in \code{rownames(sce)}. It is ignored
#' when \code{foreground = "solid"} or \code{foreground = "none"}. Defaults to
#' \code{"logcounts"}.
#'
#' @param fg_type A character value indicating how the foreground variable should
#' be represented when \code{foreground} is not \code{"solid"} or
#' \code{"none"}. Accepted values are \code{"auto"}, \code{"continuous"},
#' \code{"categorical"}, \code{"bivariate"} and \code{"combination"}. With
#' \code{"auto"}, feature values are treated as continuous and \code{colData}
#' variables are treated as continuous only when they are numeric.
#' \code{"bivariate"} blends colours from exactly two numeric foreground
#' variables.
#' \code{"combination"} calls 2 to 5 numeric foreground variables expressed
#' using \code{fg_threshold} and treats the resulting combinations as
#' categorical values. Defaults to \code{"auto"}.
#'
#' @param fg_threshold Thresholds used only when
#' \code{fg_type = "combination"}. It can be a single numeric value, a numeric
#' vector with one threshold per foreground variable, a named numeric vector
#' containing all foreground variable names, or \code{"median"}. Values are
#' expressed when \code{value > threshold}; with the default \code{0}, exact zero
#' values are not expressed. Defaults to \code{0}.
#'
#' @param fg_order A character value controlling the drawing order of foreground
#' cells within each panel. Accepted values are \code{"auto"}, \code{"input"},
#' \code{"random"}, \code{"ascending"} and \code{"descending"}.
#' \code{"auto"} uses \code{"ascending"} for bivariate and combination
#' foregrounds and \code{"input"} for other foreground types. For bivariate
#' foregrounds, ascending and descending order by the sum of the two scaled
#' foreground values. For combination foregrounds, ascending and descending
#' order by deterministic combination level, from all-negative to all-positive.
#' See \code{\link{plotOverlay}} for details. Defaults to \code{"auto"}.
#'
#' @param bg_point_size Numeric value with the point size of the background
#' layer in each panel. Defaults to \code{2}.
#'
#' @param fg_point_size Numeric value with the point size of the foreground
#' layer in each panel. Defaults to \code{0.5}.
#'
#' @param bg_palette Palette specification for the background layer. It can be
#' \code{NULL}, a palette name, a vector of colours or a palette function. When
#' \code{background = "solid"}, \code{NULL} uses \code{"gray70"} and a single
#' valid colour can be used to set the solid background colour. It is ignored
#' when \code{background = "none"}. For variable backgrounds, the value is
#' passed to \code{\link{getPalette}}. Defaults to \code{NULL}, which uses
#' \code{"gray_red"} for continuous backgrounds and \code{"scOverlay"} for
#' categorical backgrounds.
#'
#' @param fg_palette Palette specification for the foreground layer. When
#' \code{foreground = "solid"}, \code{NULL} uses \code{"red"} and a single valid
#' colour can be used to set the solid foreground colour. It is ignored when
#' \code{foreground = "none"}. For variable foregrounds, it can be \code{NULL},
#' a palette name, a vector of colours or a palette function. The value is
#' passed to \code{\link{getPalette}}. For continuous foregrounds, colour
#' vectors are interpreted as gradient anchors and interpolated internally.
#' For bivariate foregrounds, \code{NULL} uses \code{"bivariateBright"} and
#' users can provide a built-in bivariate palette name or a four-colour vector
#' with corners named \code{low_xy}, \code{high_x}, \code{high_y} and
#' \code{high_xy}. For combination foregrounds, \code{NULL} uses
#' \code{"combinationBright"}; user-supplied palettes still override this
#' default. Defaults to \code{NULL}, which uses \code{"gray_red"} for continuous
#' foregrounds and \code{"scOverlay"} for ordinary categorical foregrounds.
#'
#' @param bg_raster Logical value. If \code{TRUE}, rasterize the background point
#' layer in each panel using \code{ggrastr::geom_point_rast()}. Defaults to
#' \code{FALSE}.
#'
#' @param fg_raster Logical value. If \code{TRUE}, rasterize the foreground point
#' layer in each panel using \code{ggrastr::geom_point_rast()}. Defaults to
#' \code{FALSE}.
#'
#' @param raster_dpi Positive numeric value with the resolution, in dots per
#' inch, used for rasterized point layers. Defaults to \code{300}.
#'
#' @param bg_dimming Numeric value between 0 and 1 controlling the opacity of the
#' white dimming layer placed between background and foreground in each panel. A
#' value of \code{0} does not dim the background and a value of \code{1}
#' completely covers it with white. Defaults to \code{0.6}.
#'
#' @param fg_limits Limits for the foreground colour scale when the foreground is
#' continuous. It can be \code{NULL}, a strictly increasing numeric vector of
#' length two, or \code{"shared"}. For \code{fg_type = "bivariate"},
#' \code{NULL} uses automatic per-panel scaling and \code{"shared"} computes
#' separate shared ranges for the two numeric foreground variables. Numeric
#' foreground limits are not supported for bivariate colours. For
#' \code{fg_type = "combination"}, \code{NULL} and \code{"shared"} are accepted
#' but numeric limits are not supported. Defaults to \code{NULL}.
#'
#' @param bg_limits Limits for the background colour scale when the background is
#' continuous. It can be \code{NULL}, a strictly increasing numeric vector of
#' length two, or \code{"shared"}. With \code{"shared"}, one range is computed
#' over all cells and reused in every panel. Defaults to \code{NULL}.
#'
#' @param bg_legend Logical value indicating whether to show the background
#' legend in each panel. It is ignored when \code{background = "solid"}. Defaults
#' to \code{TRUE}.
#'
#' @param fg_legend Logical value indicating whether to show the foreground
#' legend in each panel. It is ignored when \code{foreground = "solid"}. Defaults
#' to \code{TRUE}. Bivariate foregrounds do not use a normal ggplot foreground
#' legend; use \code{\link{plotBivariateLegend}} to create a separate bivariate
#' legend.
#'
#' @param bg_legend_title Character value with the title of the background
#' legend. If \code{NULL}, \code{background} is used. Defaults to \code{NULL}.
#'
#' @param fg_legend_title Character value with the title of the foreground
#' legend. If \code{NULL}, \code{foreground} is used. For combination plots, two
#' foreground variables are joined with \code{" / "}, while larger combinations
#' use \code{"Gene combination"}. Defaults to \code{NULL}.
#'
#' @param shared_legend Logical value indicating whether compatible legends
#' should be collected into a shared patchwork legend. Defaults to \code{FALSE}.
#'
#' @return A \pkg{patchwork} object combining one \code{\link{plotOverlay}} plot
#' per group. If no valid groups produce a panel, the function returns
#' \code{NULL} with a warning.
#'
#' @seealso \code{\link{plotOverlay}},
#' \code{\link{plotOverlayPerNestedGroup}},
#' \code{\link{plotGeneOverlay}}, \code{\link{getPalette}} and
#' \code{\link{listPalettes}}.
#'
#' @examples
#' data("sce_overlay_example")
#'
#' # Plot SOX10 expression separately for each sample type.
#' p1 <- plotOverlayPerGroup(
#'     sce = sce_overlay_example,
#'     foreground = "SOX10",
#'     group_col = "sample_type",
#'     reduced_dim = "TSNE",
#'     background = "solid",
#'     fg_order = "ascending"
#' )
#' p1
#'
#' # Use the same grouping but colour the background by cluster.
#' p2 <- plotOverlayPerGroup(
#'     sce = sce_overlay_example,
#'     foreground = "S100B",
#'     group_col = "sample_type",
#'     reduced_dim = "TSNE",
#'     background = "cluster",
#'     fg_order = "ascending"
#' )
#' p2
#'
#' # Specify an explicit group order.
#' p3 <- plotOverlayPerGroup(
#'     sce = sce_overlay_example,
#'     foreground = "SOX10",
#'     group_col = "sample_type",
#'     group_order = c("PNF", "ANNUBP", "MPNST"),
#'     reduced_dim = "TSNE",
#'     background = "solid"
#' )
#' p3
#'
#' # Compare expression across samples with shared foreground limits and a
#' # collected legend.
#' p4 <- plotOverlayPerGroup(
#'     sce = sce_overlay_example,
#'     group_col = "sample",
#'     foreground = "Gene1",
#'     background = "cluster",
#'     fg_limits = "shared",
#'     shared_legend = TRUE
#' )
#' p4
#'
#' # Use shared bivariate scaling across sample panels.
#' p5 <- plotOverlayPerGroup(
#'     sce = sce_overlay_example,
#'     group_col = "sample",
#'     foreground = c("Gene1", "Gene2"),
#'     fg_type = "bivariate",
#'     fg_limits = "shared"
#' )
#' p5
#'
#' @export

plotOverlayPerGroup <- function(
		sce,
		group_col,
		foreground,
		background = "solid",
			reduced_dim = "TSNE",
			group_order = NULL,
			fg_subset_cells = NULL,
			bg_type = c("auto","continuous","categorical"),
			bg_assay = "logcounts",
			fg_assay = "logcounts",
			fg_type = c("auto","continuous","categorical","bivariate","combination"),
			fg_threshold = 0,
		fg_order = c("auto", "input", "random", "ascending", "descending"),
		bg_point_size = 2,
		fg_point_size = 0.5,
		bg_palette = NULL,
		fg_palette = NULL,
		bg_raster = FALSE,
		fg_raster = FALSE,
		raster_dpi = 300,
		bg_dimming = 0.6,
		fg_limits = NULL,
		bg_limits = NULL,
		bg_legend = TRUE,
		fg_legend = TRUE,
		bg_legend_title = NULL,
		fg_legend_title = NULL,
		shared_legend = FALSE
) {
	fg_legend_missing <- missing(fg_legend)
	.validate_sce(sce)
	.validate_coldata_column(sce, group_col, "group_col")
	.validate_raster_args(bg_raster, fg_raster, raster_dpi)
	.validate_layer_limits(bg_limits, "bg_limits")
	if (!.is_flag(shared_legend)) {
		stop("`shared_legend` must be TRUE or FALSE.", call. = FALSE)
	}
	bg_type <- match.arg(bg_type)
	fg_type <- match.arg(fg_type)
	fg_order <- match.arg(fg_order)
	if (identical(fg_type, "bivariate")) {
		.validate_bivariate_fg_limits(fg_limits)
		.handle_bivariate_fg_legend(fg_legend, fg_legend_missing)
		fg_legend <- FALSE
	} else if (identical(fg_type, "combination")) {
		.validate_combination_fg_limits(fg_limits)
	} else {
		.validate_layer_limits(fg_limits, "fg_limits")
	}
	
	vals <- colData(sce)[[group_col]]
	vals_chr <- tryCatch(as.character(vals), error = function(e) rep(NA_character_, ncol(sce)))
	if (length(vals_chr) != ncol(sce)) {
		vals_chr <- rep(NA_character_, ncol(sce))
	}
	
	if (!is.null(group_order)) {
		groups <- as.character(group_order)
	} else if (is.factor(vals)) {
		groups <- as.character(levels(vals))
	} else {
		groups <- unique(vals_chr)
	}
	
	groups <- groups[!is.na(groups) & nzchar(groups)]
	groups <- groups[!duplicated(groups)]
	
	if (length(groups) == 0) {
		warning("scOverlay - plotOverlayPerGroup: No valid groups found in '", group_col, "'. Returning NULL.")
		return(NULL)
	}
	
	subset_idx <- .eval_subset(fg_subset_cells, sce)
	layer_settings <- .resolve_grouped_layer_settings(
		sce = sce,
		background = background,
		foreground = foreground,
		bg_type = bg_type,
		fg_type = fg_type,
		bg_assay = bg_assay,
		fg_assay = fg_assay,
		fg_threshold = fg_threshold,
		subset_idx = subset_idx,
		bg_limits = bg_limits,
		fg_limits = fg_limits
	)
	
	plots <- list()
	for (g in groups) {
		idx <- vals_chr == g
		idx[is.na(idx)] <- FALSE
		panel_idx <- idx & subset_idx
		
		p <- .with_empty_subset_warning_muffled(plotOverlay(
			sce = sce,
			foreground = foreground,
			background = background,
				reduced_dim = reduced_dim,
				fg_subset_cells = panel_idx,
				bg_type = bg_type,
				bg_assay = bg_assay,
				fg_assay = fg_assay,
			fg_threshold = fg_threshold,
			fg_type = fg_type,
			fg_order = fg_order,
			bg_point_size = bg_point_size,
			fg_point_size = fg_point_size,
			bg_palette = bg_palette,
			fg_palette = fg_palette,
			bg_raster = bg_raster,
			fg_raster = fg_raster,
			raster_dpi = raster_dpi,
			bg_dimming = bg_dimming,
			fg_limits = layer_settings$fg_limits,
			bg_limits = layer_settings$bg_limits,
			bg_legend = bg_legend,
			fg_legend = fg_legend,
			bg_legend_title = bg_legend_title,
			fg_legend_title = fg_legend_title,
			title = g,
			bg_values = layer_settings$bg_values,
			fg_values = layer_settings$fg_values
		))
		plots[[g]] <- p
	}
	
	if (length(plots) == 0) {
		warning("scOverlay - plotOverlayPerGroup: Expected ", length(groups), " plots but none found. Returning NULL.")
		return(NULL)
	}
	grid_plots <- patchwork::wrap_plots(plots, ncol = length(plots))
	if (isTRUE(shared_legend)) {
		grid_plots <- grid_plots + patchwork::plot_layout(guides = "collect")
	}
	return(grid_plots)
}
