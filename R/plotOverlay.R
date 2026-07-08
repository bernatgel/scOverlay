#' Plot layered information on a single-cell embedding
#'
#' \code{plotOverlay} creates a layered plot of a reduced dimension
#' representation stored in a \linkS4class{SingleCellExperiment} object. It
#' first creates a background layer, dims this background with a semi-transparent
#' white overlay, and finally draws a foreground variable on top using an
#' independent colour scale. When a layer is set to \code{"none"}, the
#' corresponding point layer is kept in the plot structure without drawing
#' visible points.
#'
#' This representation is useful to highlight gene expression, scores,
#' annotations or selected groups of cells while keeping the complete embedding
#' visible as context.
#'
#' @details
#' The function uses two visual layers. The background layer is used to show the
#' global structure of the embedding, usually coloured by a cell annotation such
#' as a cluster, sample or cell type. Either layer can also be set to
#' \code{"none"} to keep the layer position without drawing visible points. The
#' foreground layer is then drawn on top and can represent a fixed-colour
#' selection, a feature from one of the assays, typically a gene expression
#' value, or a variable stored in \code{colData(sce)}.
#'
#' The foreground can be restricted to a subset of cells using
#' \code{fg_subset_cells}. This only affects the foreground layer: cells not
#' selected by \code{fg_subset_cells} are still shown in the background. This is
#' the main idea of the overlay plot and allows highlighting a particular group
#' of cells without losing the context of the whole dataset.
#'
#' The foreground colour scale is independent from the background colour scale.
#' Palettes can be specified by name, as a vector of colours, or as a palette
#' function. See \code{\link{getPalette}} and \code{\link{listPalettes}} for the
#' available palettes and palette handling rules. For continuous palettes, colour
#' vectors are used as anchors and interpolated internally.
#'
#' For dense embeddings, the background and foreground point layers can be
#' rasterized independently with \code{bg_raster} and \code{fg_raster}. Only the
#' point layers are rasterized; axes, labels, legends and other plot elements
#' remain vector graphics.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object containing the reduced
#' dimensions and the data to plot. The object must contain the reduced
#' dimension specified in \code{reduced_dim}. If \code{background} or
#' \code{foreground} is a feature, the object must also contain the
#' corresponding assay specified in \code{bg_assay} or \code{fg_assay}.
#'
#' @param foreground A character value specifying the variable to draw in the
#' foreground layer. It can be \code{"solid"} to draw selected foreground cells
#' with a fixed colour, \code{"none"} to draw no visible foreground points, a
#' feature in \code{rownames(sce)}, in which case values are taken from
#' \code{fg_assay}, or a column in \code{colData(sce)}. When
#' \code{fg_type = "bivariate"}, it must be a character vector of exactly two
#' numeric features or \code{colData(sce)} columns. When
#' \code{fg_type = "combination"}, it must contain between 2 and 5 numeric
#' features or \code{colData(sce)} columns. This argument has no default and
#' must be provided by the user.
#'
#' @param background A single character value specifying the background layer.
#' It can be \code{"solid"} to draw all background cells with a fixed colour,
#' \code{"none"} to draw no visible background points, a feature in
#' \code{rownames(sce)}, or a column in \code{colData(sce)}. Defaults to
#' \code{"solid"}.
#'
#' @param reduced_dim A single character value with the name of the reduced
#' dimension to plot. It must be one of \code{reducedDimNames(sce)}. Only the
#' first two dimensions are used. Defaults to \code{"TSNE"}.
#'
#' @param fg_subset_cells A cell subset defining which cells are shown in the
#' foreground layer. It can be \code{NULL}, a logical vector of length
#' \code{ncol(sce)}, a character vector with cell names, a function receiving
#' \code{sce} and returning a logical vector, or an expression evaluated in
#' \code{colData(sce)}. \code{NULL} selects all cells. \code{NA} values are
#' treated as \code{FALSE}. If no cells match, an empty foreground layer is
#' drawn and a warning is emitted. Defaults to \code{NULL}.
#'
#' @param bg_type A character value indicating how the background variable should
#' be represented when \code{background} is not \code{"solid"} or \code{"none"}.
#' Accepted values are \code{"auto"}, \code{"continuous"} and
#' \code{"categorical"}. With \code{"auto"}, feature values are treated as
#' continuous and \code{colData} variables are treated as continuous only when
#' they are numeric. Defaults to \code{"auto"}.
#'
#' @param fg_type A character value indicating how the foreground variable should
#' be represented when \code{foreground} is not \code{"solid"} or \code{"none"}.
#' Accepted values are \code{"auto"}, \code{"continuous"},
#' \code{"categorical"}, \code{"bivariate"} and \code{"combination"}. With
#' \code{"auto"}, feature values are treated as continuous and \code{colData}
#' variables are treated as continuous only when they are numeric.
#' \code{"bivariate"} blends colours from two continuous foreground variables.
#' \code{"combination"} calls each of 2 to 5 numeric foreground variables
#' expressed when its value is greater than \code{fg_threshold}, then treats the
#' expressed/not-expressed combinations as categorical values. Defaults to
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
#' @param fg_threshold Thresholds used only when
#' \code{fg_type = "combination"}. It can be a single numeric value, a numeric
#' vector with one threshold per foreground variable, a named numeric vector
#' containing all foreground variable names, or the single character value
#' \code{"median"}. Values are considered expressed when
#' \code{value > threshold}; with the default \code{0}, exact zero values are
#' not expressed. Defaults to \code{0}.
#'
#' @param fg_order A character value controlling the drawing order of foreground
#' cells. Accepted values are \code{"input"}, \code{"random"},
#' \code{"ascending"} and \code{"descending"}. \code{"input"} keeps the current
#' order, \code{"random"} randomly shuffles foreground cells,
#' \code{"ascending"} draws lower foreground values first and higher values on
#' top, and \code{"descending"} does the opposite. \code{"ascending"} is useful
#' for expression plots because cells with higher expression are drawn last. For
#' reproducible random ordering, call \code{set.seed()} before
#' \code{plotOverlay()}. Defaults to \code{"input"}.
#'
#' @param bg_point_size Numeric value with the point size of the background
#' layer. Defaults to \code{4}.
#'
#' @param fg_point_size Numeric value with the point size of the foreground
#' layer. Defaults to \code{1}.
#'
#' @param bg_palette Palette specification for the background layer. It is
#' ignored when \code{background = "none"}. When \code{background = "solid"},
#' \code{NULL} uses \code{"gray70"} and a single valid colour can be used to set
#' the solid background colour. For variable backgrounds, it can be \code{NULL},
#' a palette name, a vector of colours or a palette function, and the value is
#' passed to \code{\link{getPalette}}. Defaults to \code{NULL}, which uses
#' \code{"gray_red"} for continuous backgrounds and \code{"scOverlay"} for
#' categorical backgrounds.
#'
#' @param fg_palette Palette specification for the foreground layer. It is
#' ignored when \code{foreground = "none"}. When \code{foreground = "solid"},
#' \code{NULL} uses \code{"red"} and a single valid colour can be used to set the
#' solid foreground colour. For variable foregrounds, it can be \code{NULL}, a
#' palette name, a vector of colours or a palette function, and the value is
#' passed to \code{\link{getPalette}}. For continuous foregrounds, colour vectors
#' are interpreted as gradient anchors and interpolated internally. For
#' bivariate foregrounds, \code{NULL} uses \code{"bivariateBright"} and users
#' can provide a built-in bivariate palette name or a four-colour vector with
#' corners named \code{low_xy}, \code{high_x}, \code{high_y} and
#' \code{high_xy}. The first foreground variable maps to the x direction and
#' \code{high_x}; the second maps to the y direction and \code{high_y}. For
#' combination foregrounds, \code{NULL} uses \code{"combinationBright"};
#' user-supplied palettes still override this default. Defaults to \code{NULL},
#' which uses \code{"gray_red"} for continuous foregrounds and \code{"scOverlay"}
#' for ordinary categorical foregrounds.
#'
#' @param bg_raster Logical value. If \code{TRUE}, rasterize the background point
#' layer using \code{ggrastr::geom_point_rast()}. Defaults to \code{FALSE}.
#'
#' @param fg_raster Logical value. If \code{TRUE}, rasterize the foreground point
#' layer using \code{ggrastr::geom_point_rast()}. Defaults to \code{FALSE}.
#'
#' @param raster_dpi Positive numeric value with the resolution, in dots per
#' inch, used for rasterized point layers. Defaults to \code{300}.
#'
#' @param bg_dimming Numeric value between 0 and 1 controlling the opacity of the
#' white dimming layer placed between background and foreground. A value of
#' \code{0} does not dim the background and a value of \code{1} completely covers
#' it with white. Defaults to \code{0.9}.
#'
#' @param fg_limits Limits for the foreground colour scale when the foreground is
#' continuous. It can be \code{NULL}, a strictly increasing numeric vector of
#' length two, or \code{"shared"}. For \code{fg_type = "bivariate"},
#' \code{NULL} uses automatic scaling and \code{"shared"} can be used by grouped
#' plotting helpers for comparable bivariate colours across panels. Numeric
#' foreground limits are not supported for bivariate colours. For
#' \code{fg_type = "combination"}, \code{NULL} and \code{"shared"} are accepted
#' but numeric limits are not supported. It is ignored with a warning when the
#' foreground layer is not continuous, bivariate or combination. Defaults to
#' \code{NULL}.
#'
#' @param bg_limits Limits for the background colour scale when the background is
#' continuous. It can be \code{NULL}, a strictly increasing numeric vector of
#' length two, or \code{"shared"}. In a single plot, \code{"shared"} resolves to
#' the range of background values in that plot. It is ignored with a warning when
#' the background layer is not continuous. Defaults to \code{NULL}.
#'
#' @param bg_legend Logical value indicating whether to show the background
#' legend. It is ignored when \code{background = "solid"} or
#' \code{background = "none"}. Defaults to \code{TRUE}.
#'
#' @param fg_legend Logical value indicating whether to show the foreground
#' legend. It is ignored when \code{foreground = "solid"} or
#' \code{foreground = "none"}. Bivariate foregrounds do not use a normal ggplot
#' foreground legend; use \code{\link{plotBivariateLegend}} to create a separate
#' bivariate legend. Defaults to \code{TRUE}.
#'
#' @param bg_legend_title Character value with the title of the background
#' legend. If \code{NULL}, \code{background} is used. Defaults to \code{NULL}.
#'
#' @param fg_legend_title Character value with the title of the foreground
#' legend. If \code{NULL}, \code{foreground} is used. For combination plots, two
#' foreground variables are joined with \code{" / "}, while larger combinations
#' use \code{"Gene combination"}. Defaults to \code{NULL}.
#'
#' @param title Character value with the plot title. If \code{NULL},
#' \code{foreground} is used, except when \code{foreground = "none"}, where an
#' empty title is used. Defaults to \code{NULL}.
#'
#' @param bg_values Categorical value set for the background
#' scale. This is used by grouped plotting helpers to keep categorical colours
#' consistent across panels. Not intended for users.
#'
#' @param fg_values Categorical value set for the foreground
#' scale. This is used by grouped plotting helpers to keep categorical colours
#' consistent across panels. Not intended for users.
#'
#' @return A \pkg{ggplot2} object with the background, dimming and foreground
#' layers.
#'
#' @seealso \code{\link{plotGeneOverlay}},
#' \code{\link{plotOverlayPerGroup}},
#' \code{\link{plotOverlayPerNestedGroup}}, \code{\link{getPalette}} and
#' \code{\link{listPalettes}}.
#'
#' @examples
#' data("sce_overlay_example")
#'
#' # Overlay the expression of a gene on top of the cluster annotation.
#' p1 <- plotOverlay(
#'     sce = sce_overlay_example,
#'     reduced_dim = "TSNE",
#'     background = "cluster",
#'     foreground = "SOX10"
#' )
#' p1
#'
#' # Use a solid background and draw higher expression values on top.
#' p2 <- plotOverlay(
#'     sce = sce_overlay_example,
#'     reduced_dim = "TSNE",
#'     background = "solid",
#'     foreground = "S100B",
#'     fg_order = "ascending"
#' )
#' p2
#'
#' # Restrict the foreground to one sample type while keeping all cells in the
#' # background.
#' p3 <- plotOverlay(
#'     sce = sce_overlay_example,
#'     reduced_dim = "UMAP",
#'     background = "cluster",
#'     foreground = "SOX10",
#'     fg_subset_cells = quote(sample_type == "ANNUBP")
#' )
#' p3
#'
#' # Use a custom continuous palette and rasterize the background point layer.
#' p4 <- plotOverlay(
#'     sce = sce_overlay_example,
#'     reduced_dim = "TSNE",
#'     background = "sample_type",
#'     foreground = "SOX10",
#'     fg_palette = c("#E6E6E6", "#EF6F6C", "#67000D"),
#'     bg_raster = TRUE
#' )
#' p4
#'
#' @export

plotOverlay <- function(
		sce,
		foreground,                    
		background = "solid",          
		reduced_dim = "TSNE",          
		fg_subset_cells = NULL,        
		bg_type = c("auto","continuous","categorical"),
		fg_type = c("auto","continuous","categorical","bivariate","combination"),
		bg_assay = "logcounts",        
		fg_assay = "logcounts",        
		fg_threshold = 0,
		fg_order = c("input", "random", "ascending", "descending"),
		bg_point_size = 4,
		fg_point_size = 1,
		bg_palette = NULL,
		fg_palette = NULL,
		bg_raster = FALSE,
		fg_raster = FALSE,
		raster_dpi = 300,
		bg_dimming = 0.9,              
		fg_limits = NULL,
		bg_limits = NULL,
		bg_legend = TRUE,
		fg_legend = TRUE,
		bg_legend_title = NULL,
		fg_legend_title = NULL,
		title = NULL,
		bg_values = NULL,
		fg_values = NULL
) {
	fg_legend_missing <- missing(fg_legend)
	.validate_sce(sce)
	if (!is.numeric(bg_dimming) || length(bg_dimming) != 1L ||
			is.na(bg_dimming) || !is.finite(bg_dimming) ||
			bg_dimming < 0 || bg_dimming > 1) {
		stop("`bg_dimming` must be a single numeric value between 0 and 1.", call. = FALSE)
	}
	if (!.is_flag(bg_legend)) {
		stop("`bg_legend` must be TRUE or FALSE.", call. = FALSE)
	}
	if (!.is_flag(fg_legend)) {
		stop("`fg_legend` must be TRUE or FALSE.", call. = FALSE)
	}
	.validate_raster_args(bg_raster, fg_raster, raster_dpi)
	.validate_string_arg(background, "background")
	.validate_layer_limits(bg_limits, "bg_limits")
	
	# Extract embedding coordinates
	coords <- .get_reduced_dim_coords(sce, reduced_dim)
	
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

		bg <- .resolve_overlay_variable(
			sce = sce,
			variable = background,
			layer_name = "background",
			assay = bg_assay,
			type = bg_type
		)
		coords$background <- bg$values
		if (identical(bg$type, "categorical")) {
			coords$background <- as.character(coords$background)
		}
		
		if (identical(fg_type, "bivariate")) {
			fg <- .prepare_bivariate_values(
				sce = sce,
				features = foreground,
				assay = fg_assay,
				limits = fg_limits,
				palette = .resolve_fg_palette(fg_palette, "bivariate")
			)
		} else if (identical(fg_type, "combination")) {
			fg <- .prepare_combination_values(
				sce = sce,
				features = foreground,
				assay = fg_assay,
				threshold = fg_threshold
			)
		} else {
			.validate_string_arg(foreground, "foreground")
			fg <- .resolve_overlay_variable(
			sce = sce,
			variable = foreground,
					layer_name = "foreground",
					assay = fg_assay,
					type = fg_type
				)
		}
		coords$foreground <- if (identical(fg$type, "bivariate")) {
			fg$colours
		} else if (identical(fg$type, "combination")) {
			fg$labels
		} else {
			fg$values
		}
		if (identical(fg$type, "bivariate")) {
			coords$.fg_x <- fg$values_x
			coords$.fg_y <- fg$values_y
			coords$.fg_colour <- fg$colours
		}
		fg_type <- fg$type
		
	if (identical(fg$type, "none")) {
		coords_fg <- coords[0, , drop = FALSE]
	} else {
		# Apply subset filter to foreground only (background remains complete)
		idx <- .eval_subset(fg_subset_cells, sce)
		if (identical(fg$type, "bivariate")) {
			coords$foreground[!idx] <- NA_character_
			coords$.fg_colour[!idx] <- NA_character_
		} else {
			coords$foreground[!idx] <- NA
		}

		# Create filtered coords for foreground layer (only subset cells with non-NA foreground)
		coords_fg <- coords[!is.na(coords$foreground), ]
		coords_fg <- .order_foreground(coords_fg, fg_order)
	}

	bg_limits <- .resolve_layer_limits(
		limits = bg_limits,
		values = coords$background,
		type = bg$type,
		arg_name = "bg_limits",
		layer_name = "background"
	)
	if (identical(fg$type, "bivariate")) {
		fg_limits <- list(x = fg$limits_x, y = fg$limits_y)
	} else if (identical(fg$type, "combination")) {
		fg_limits <- NULL
	} else {
		fg_limits <- .resolve_layer_limits(
			limits = fg_limits,
			values = coords_fg$foreground,
			type = fg$type,
			arg_name = "fg_limits",
			layer_name = "foreground"
		)
	}
	
		# === Build ggplot ===
		# Background layer
		if (identical(bg$type, "none")) {
			p <- ggplot(coords, aes(.data$X, .data$Y)) +
				.geom_point_overlay(
					data = coords[0, , drop = FALSE],
					mapping = aes(.data$X, .data$Y),
					raster = FALSE,
					color = "transparent",
					size = 0,
					na.rm = TRUE
				)
		} else if (identical(bg$type, "solid")) {
			bg_colour <- .resolve_solid_bg_colour(bg_palette)
			p <- ggplot(coords, aes(.data$X, .data$Y)) +
				.geom_point_overlay(
				raster = bg_raster,
				raster_dpi = raster_dpi,
				color = bg_colour,
					size = bg_point_size
				)
		} else if (identical(bg$type, "continuous")) {
			pal_bg <- getPalette(bg_palette, type = "continuous", n = 256)
			p <- ggplot(coords, aes(.data$X, .data$Y, colour = .data$background)) +
				.geom_point_overlay(
					raster = bg_raster,
					raster_dpi = raster_dpi,
					size = bg_point_size
				) +
				scale_colour_gradientn(
					name = bg_legend_title %||% background,
					colours = pal_bg,
					limits = bg_limits,
					guide = if (bg_legend) "colourbar" else "none"
				)
		} else {
			bg_values <- .resolve_categorical_values(bg_values, coords$background)
			pal_bg <- getPalette(
				bg_palette,
			type = "categorical",
			n = max(1L, length(bg_values)),
			values = bg_values
		)
		p <- ggplot(coords, aes(.data$X, .data$Y, colour = .data$background)) +
			.geom_point_overlay(
				raster = bg_raster,
				raster_dpi = raster_dpi,
				size = bg_point_size
			) +
			scale_colour_manual(
				name = bg_legend_title %||% background,
				values = pal_bg,
				limits = bg_values,
				guide = if (bg_legend) "legend" else "none"
			)
	}
	
	# Dim background
	dim_alpha <- if (identical(bg$type, "none")) 0 else bg_dimming
	p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
										fill = "white", alpha = dim_alpha)
	
		# Foreground layer
			if (fg$type == "none") {
				p <- p + .geom_point_overlay(
					data = coords_fg,
					mapping = aes(.data$X, .data$Y),
					raster = FALSE,
					color = "transparent",
					size = 0,
					na.rm = TRUE
				)
			} else if (fg$type == "solid") {
				fg_colour <- .resolve_solid_fg_colour(fg_palette)
				p <- p + .geom_point_overlay(
					data = coords_fg,
					mapping = aes(.data$X, .data$Y),
					raster = fg_raster,
					raster_dpi = raster_dpi,
					color = fg_colour,
					size = fg_point_size,
					na.rm = TRUE
				)
			} else {
				p <- p + ggnewscale::new_scale_color()
			}
		
			if (fg$type == "continuous") {
			pal_fg <- getPalette(fg_palette, type = "continuous", n = 256)

		p <- p + .geom_point_overlay(
				data = coords_fg,
				mapping = aes(.data$X, .data$Y, colour = .data$foreground),
				raster = fg_raster,
				raster_dpi = raster_dpi,
				size = fg_point_size,
				na.rm = TRUE
			) +
			scale_colour_gradientn(
				name = fg_legend_title %||% foreground,
				colours = pal_fg,
				limits = fg_limits,
				guide = if (fg_legend) "colourbar" else "none"
			)
		
			} else if (fg$type == "categorical" || fg$type == "combination") {
			if (identical(fg$type, "combination") && is.null(fg_values)) {
				fg_values <- levels(fg$labels)
			}
			coords_fg$foreground <- as.character(coords_fg$foreground)
			fg_values <- .resolve_categorical_values(fg_values, coords_fg$foreground)
			n_fg <- max(1L, length(fg_values))
			pal_fg <- getPalette(
			.resolve_fg_palette(fg_palette, fg$type),
			type = "categorical",
			n = n_fg,
			values = fg_values
		)

		p <- p + .geom_point_overlay(
				data = coords_fg,
				mapping = aes(.data$X, .data$Y, colour = .data$foreground),
				raster = fg_raster,
				raster_dpi = raster_dpi,
				size = fg_point_size,
				na.rm = TRUE
			) +
			scale_colour_manual(
				name = fg_legend_title %||% .foreground_legend_title(foreground, fg$type),
				values = pal_fg,
				limits = fg_values,
				guide = if (fg_legend) "legend" else "none"
			)
	} else if (fg$type == "bivariate") {
		# Bivariate colours are precomputed from two continuous variables, so no
		# ggplot scale or automatic legend is added in this first implementation.
		p <- p + .geom_point_overlay(
			data = coords_fg,
			mapping = aes(.data$X, .data$Y),
			raster = fg_raster,
			raster_dpi = raster_dpi,
			colour = coords_fg$.fg_colour,
			size = fg_point_size,
			na.rm = TRUE
		)
	}
	
	# Titles
	p <- p + theme_classic()
	if (identical(bg$type, "none")) {
		p <- p + ggplot2::coord_cartesian(
			xlim = range(coords$X, na.rm = TRUE),
			ylim = range(coords$Y, na.rm = TRUE)
		)
	}
	if (!is.null(title)) p <- p + labs(title = title)
	else if (identical(fg$type, "none")) p <- p + labs(title = "")
	else p <- p + labs(title = paste(fg$label, collapse = " + "))
	
	return(p)
}

#' Plot a standalone bivariate colour legend
#'
#' \code{plotBivariateLegend} creates a square tile legend for bivariate
#' foreground palettes used by \code{\link{plotOverlay}} when
#' \code{fg_type = "bivariate"}. It returns a standalone \pkg{ggplot2} object
#' that can be combined manually with an overlay plot, for example with
#' \pkg{patchwork}.
#'
#' @param features Character vector of length two. The first feature labels the
#' x direction and maps to the \code{high_x} palette corner; the second feature
#' labels the y direction and maps to \code{high_y}.
#'
#' @param palette Bivariate palette specification. It can be one of the built-in
#' bivariate palette names, such as \code{"bivariateBright"}, or a four-colour
#' vector with corners named \code{low_xy}, \code{high_x}, \code{high_y} and
#' \code{high_xy}. Unnamed four-colour vectors are interpreted in that order.
#' Defaults to \code{"bivariateBright"}.
#'
#' @param n_bins Number of bins per axis in the displayed colour grid. Must be a
#' single integer greater than or equal to 2. Defaults to \code{9}.
#'
#' @return A \pkg{ggplot2} object containing the bivariate colour grid.
#'
#' @examples
#' leg <- plotBivariateLegend(c("SOX10", "MKI67"))
#' leg
#'
#' @export

plotBivariateLegend <- function(features,
																palette = "bivariateBright",
																n_bins = 9L) {
	return(.plot_bivariate_palette_grid(
		palette = palette,
		features = features,
		n_bins = n_bins,
		title = NULL
	))
}

.plot_bivariate_palette_grid <- function(palette,
																				 features = c("Feature 1", "Feature 2"),
																				 n_bins = 9L,
																				 title = NULL) {
	if (!is.character(features) || length(features) != 2L || anyNA(features)) {
		stop("`features` must be a character vector of length 2.", call. = FALSE)
	}
	if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) || n_bins < 2L) {
		stop("`n_bins` must be a single integer greater than or equal to 2.", call. = FALSE)
	}
	n_bins <- as.integer(n_bins)
	palette_matrix <- .bivariate_palette(palette, n_bins)
	legend_data <- expand.grid(
		x = seq_len(n_bins),
		y = seq_len(n_bins),
		KEEP.OUT.ATTRS = FALSE
	)
	legend_data$colour <- as.vector(t(palette_matrix))
	p <- ggplot2::ggplot(legend_data, ggplot2::aes(.data$x, .data$y, fill = .data$colour)) +
		ggplot2::geom_tile() +
		ggplot2::scale_fill_identity() +
		ggplot2::scale_x_continuous(
			name = features[[1]],
			breaks = c(1, n_bins),
			labels = c("Low", "High"),
			expand = c(0, 0)
		) +
		ggplot2::scale_y_continuous(
			name = features[[2]],
			breaks = c(1, n_bins),
			labels = c("Low", "High"),
			expand = c(0, 0)
		) +
		ggplot2::labs(x = features[[1]], y = features[[2]]) +
		ggplot2::coord_fixed() +
		ggplot2::theme_classic()
	if (!is.null(title)) {
		p <- p + ggplot2::ggtitle(title)
	}
	return(p)
}



# Helper for NULL defaults
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper for scalar logical flags
.is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

# Helper to validate scalar character arguments
.validate_string_arg <- function(x, arg_name) {
	if (!is.character(x) || length(x) != 1L || is.na(x)) {
		stop("`", arg_name, "` must be a single non-NA character string.", call. = FALSE)
	}
	return(invisible(TRUE))
}

# Helper to validate colour scale limit arguments.
.validate_layer_limits <- function(limits, arg_name) {
	if (is.null(limits)) {
		return(invisible(TRUE))
	}
	if (is.character(limits) && length(limits) == 1L && !is.na(limits) &&
			identical(limits, "shared")) {
		return(invisible(TRUE))
	}
	if (!is.numeric(limits) || length(limits) != 2L) {
		stop("`", arg_name, "` must be NULL, \"shared\", or a numeric vector of length 2.", call. = FALSE)
	}
	if (any(is.na(limits)) || any(!is.finite(limits))) {
		stop("`", arg_name, "` must contain two finite non-NA values.", call. = FALSE)
	}
	if (!(limits[1] < limits[2])) {
		stop("`", arg_name, "` must be strictly increasing.", call. = FALSE)
	}
	return(invisible(TRUE))
}

# Helper to validate user-facing bivariate foreground limits.
.validate_bivariate_fg_limits <- function(limits) {
	if (is.null(limits)) {
		return(invisible(TRUE))
	}
	if (inherits(limits, "scOverlay_bivariate_limits")) {
		return(invisible(TRUE))
	}
	if (is.character(limits) && length(limits) == 1L && !is.na(limits) &&
			identical(limits, "shared")) {
		return(invisible(TRUE))
	}
	stop(
		"`fg_limits` must be NULL or \"shared\" when `fg_type = \"bivariate\"`. ",
		"Numeric foreground limits are not supported for bivariate colours.",
		call. = FALSE
	)
}

.new_bivariate_limits <- function(x, y) {
	structure(list(x = x, y = y), class = "scOverlay_bivariate_limits")
}

# Helper to choose mode-specific foreground palette defaults.
.resolve_fg_palette <- function(fg_palette, fg_type) {
	if (is.null(fg_palette) && identical(fg_type, "bivariate")) {
		return("bivariateBright")
	}
	if (is.null(fg_palette) && identical(fg_type, "combination")) {
		return("combinationBright")
	}
	return(fg_palette)
}

# Helper to warn users away from unsupported one-dimensional bivariate guides.
.handle_bivariate_fg_legend <- function(fg_legend, fg_legend_missing) {
	if (!isTRUE(fg_legend_missing) && isTRUE(fg_legend)) {
		warning(
			"`fg_legend = TRUE` is not supported for `fg_type = \"bivariate\"`. ",
			"Use `plotBivariateLegend()` to create a separate bivariate legend.",
			call. = FALSE
		)
	}
	return(invisible(FALSE))
}

# Helper to choose default foreground legend titles.
.foreground_legend_title <- function(foreground, fg_type) {
	if (identical(fg_type, "combination")) {
		if (length(foreground) == 2L) {
			return(paste(foreground, collapse = " / "))
		}
		return("Gene combination")
	}
	return(foreground)
}

# Helper to validate combination foreground limits.
.validate_combination_fg_limits <- function(limits) {
	if (is.null(limits)) {
		return(invisible(TRUE))
	}
	if (is.character(limits) && length(limits) == 1L && !is.na(limits) &&
			identical(limits, "shared")) {
		return(invisible(TRUE))
	}
	stop("`fg_limits` is not supported when `fg_type = \"combination\"`.", call. = FALSE)
}

# Helper to resolve NULL/numeric/"shared" limits for one layer.
.resolve_layer_limits <- function(limits, values, type, arg_name, layer_name, warn = TRUE) {
	.validate_layer_limits(limits, arg_name)
	if (is.null(limits)) {
		return(NULL)
	}
	if (!identical(type, "continuous")) {
		if (isTRUE(warn)) {
			warning(
				"`", arg_name, "` is ignored because the ",
				layer_name,
				" layer is not continuous.",
				call. = FALSE
			)
		}
		return(NULL)
	}
	if (is.character(limits) && identical(limits, "shared")) {
		return(.range_or_null(values, arg_name))
	}
	return(limits)
}

# Helper to compute finite ranges without breaking empty pipeline-safe plots.
.range_or_null <- function(values, arg_name) {
	values <- values[!is.na(values) & is.finite(values)]
	if (length(values) == 0L) {
		warning(
			"`", arg_name, " = \"shared\"` could not be computed because no finite values are available.",
			call. = FALSE
		)
		return(NULL)
	}
	limits <- range(values)
	if (identical(limits[1], limits[2])) {
		limits <- limits + c(-0.5, 0.5)
	}
	return(limits)
}

# Helper to resolve one numeric foreground feature or colData column.
.resolve_numeric_foreground_feature <- function(sce, feature, assay, type_label) {
	.validate_string_arg(feature, "foreground")
	if (feature %in% SummarizedExperiment::rownames(sce)) {
		.validate_string_arg(assay, "fg_assay")
		if (!(assay %in% SummarizedExperiment::assayNames(sce))) {
			stop(
				"`fg_assay` must be one of assayNames(sce). Unknown assay: ",
				assay,
				".",
				call. = FALSE
			)
		}
		return(as.numeric(SummarizedExperiment::assay(sce, assay)[feature, ]))
	}
	if (feature %in% colnames(SummarizedExperiment::colData(sce))) {
		values <- SummarizedExperiment::colData(sce)[[feature]]
		if (!is.numeric(values)) {
			if (identical(type_label, "bivariate")) {
				stop("`fg_type = 'bivariate'` requires both foreground variables to be numeric.", call. = FALSE)
			}
			stop("`fg_type = '", type_label, "'` requires foreground variables to be numeric.", call. = FALSE)
		}
		return(values)
	}
	stop(
		"`foreground` entries for `fg_type = '", type_label, "'` must be features in rownames(sce) ",
		"or numeric columns in colData(sce). Unknown foreground: ",
		feature,
		".",
		call. = FALSE
	)
}

# Helper to resolve one bivariate foreground feature or colData column.
.resolve_bivariate_feature <- function(sce, feature, assay) {
	return(.resolve_numeric_foreground_feature(sce, feature, assay, "bivariate"))
}

# Helper to validate and expand combination thresholds.
.resolve_combination_thresholds <- function(threshold, features, values) {
	if (is.character(threshold)) {
		if (length(threshold) == 1L && !is.na(threshold) && identical(threshold, "median")) {
			thresholds <- apply(values, 2L, stats::median, na.rm = TRUE)
			thresholds <- as.numeric(thresholds)
			names(thresholds) <- features
			return(thresholds)
		}
		stop("`fg_threshold` must be numeric or \"median\" when `fg_type = \"combination\"`.", call. = FALSE)
	}
	if (!is.numeric(threshold) || anyNA(threshold) || any(!is.finite(threshold))) {
		stop("`fg_threshold` must be numeric or \"median\" when `fg_type = \"combination\"`.", call. = FALSE)
	}
	if (!is.null(names(threshold)) && any(nzchar(names(threshold)))) {
		if (!all(features %in% names(threshold))) {
			stop("Named `fg_threshold` values must include all foreground variables.", call. = FALSE)
		}
		thresholds <- as.numeric(threshold[features])
		names(thresholds) <- features
		return(thresholds)
	}
	if (length(threshold) == 1L) {
		thresholds <- rep(as.numeric(threshold), length(features))
		names(thresholds) <- features
		return(thresholds)
	}
	if (length(threshold) == length(features)) {
		thresholds <- as.numeric(threshold)
		names(thresholds) <- features
		return(thresholds)
	}
	stop("`fg_threshold` must have length 1 or the same length as `foreground`.", call. = FALSE)
}

# Helper to generate all combination labels in deterministic binary order.
.combination_levels <- function(features) {
	k <- length(features)
	grid <- expand.grid(rep(list(c(FALSE, TRUE)), k), KEEP.OUT.ATTRS = FALSE)
	levels <- apply(grid, 1L, function(row) {
		paste0(features, ifelse(row, "+", "-"), collapse = "/")
	})
	return(levels)
}

# Helper to construct expression-combination foreground labels.
.prepare_combination_values <- function(sce,
																			 features,
																			 assay,
																			 threshold = 0) {
	if (!is.character(features) || anyNA(features) ||
			length(features) < 2L || length(features) > 5L) {
		stop("`fg_type = \"combination\"` supports between 2 and 5 foreground variables.", call. = FALSE)
	}
	values <- vapply(
		features,
		function(feature) .resolve_numeric_foreground_feature(sce, feature, assay, "combination"),
		numeric(ncol(sce))
	)
	if (is.null(dim(values))) {
		values <- matrix(values, ncol = length(features), dimnames = list(colnames(sce), features))
	}
	colnames(values) <- features
	thresholds <- .resolve_combination_thresholds(threshold, features, values)
	expressed <- sweep(values, 2L, thresholds, FUN = ">")
	has_na <- rowSums(is.na(values)) > 0L
	labels <- apply(expressed, 1L, function(row) {
		paste0(features, ifelse(row, "+", "-"), collapse = "/")
	})
	labels[has_na] <- NA_character_
	levels <- .combination_levels(features)
	labels <- factor(labels, levels = levels)
	return(list(
		features = features,
		values = values,
		thresholds = thresholds,
		expressed = expressed,
		labels = labels,
		type = "combination",
		source = "combination",
		label = features
	))
}

# Helper to construct bivariate foreground colours from two numeric variables.
.prepare_bivariate_values <- function(sce,
																			features,
																			assay,
																			limits = NULL,
																			n_bins = 9L,
																			palette = "bivariateBright") {
	if (!is.character(features) || length(features) != 2L || anyNA(features)) {
		stop("`fg_type = 'bivariate'` requires `foreground` to be a character vector of length 2.", call. = FALSE)
	}
	.validate_bivariate_fg_limits(limits)
	x <- .resolve_bivariate_feature(sce, features[1], assay)
	y <- .resolve_bivariate_feature(sce, features[2], assay)
	if (!is.numeric(x) || !is.numeric(y)) {
		stop("`fg_type = 'bivariate'` requires both foreground variables to be numeric.", call. = FALSE)
	}
	limits_x <- if (inherits(limits, "scOverlay_bivariate_limits")) limits$x else .range_or_null(x, "fg_limits$x")
	limits_y <- if (inherits(limits, "scOverlay_bivariate_limits")) limits$y else .range_or_null(y, "fg_limits$y")
	palette_matrix <- .bivariate_palette(palette, n_bins)
	x_scaled <- scales::rescale(x, from = limits_x)
	y_scaled <- scales::rescale(y, from = limits_y)
	x_scaled <- pmin(pmax(x_scaled, 0), 1)
	y_scaled <- pmin(pmax(y_scaled, 0), 1)
	breaks <- seq(0, 1, length.out = n_bins + 1L)
	x_bin <- cut(x_scaled, breaks = breaks, include.lowest = TRUE, labels = FALSE)
	y_bin <- cut(y_scaled, breaks = breaks, include.lowest = TRUE, labels = FALSE)
	colours <- rep(NA_character_, length(x))
	ok <- !is.na(x_bin) & !is.na(y_bin)
	colours[ok] <- palette_matrix[cbind(y_bin[ok], x_bin[ok])]
	return(list(
		feature_x = features[1],
		feature_y = features[2],
		values_x = x,
		values_y = y,
		scaled_x = x_scaled,
		scaled_y = y_scaled,
		colours = colours,
		limits_x = limits_x,
		limits_y = limits_y,
		palette = palette_matrix,
		values = colours,
		type = "bivariate",
		source = "bivariate",
		label = features
	))
}

.bivariate_palette_registry <- function() {
	return(list(
		bivariateBright = c(
			low_xy = "#e6e6e6",
			high_x = "#ff4346",
			high_y = "#55b2ff",
			high_xy = "#7a1fa2"
		),
		bivariateOrangeBlue = c(
			low_xy = "#e8e8e8",
			high_x = "#f46d43",
			high_y = "#4393c3",
			high_xy = "#542788"
		),
		bivariatePinkGreen = c(
			low_xy = "#e8e8e8",
			high_x = "#e7298a",
			high_y = "#66a61e",
			high_xy = "#762a83"
		),
		bivariateRedCyan = c(
			low_xy = "#e6e6e6",
			high_x = "#e41a1c",
			high_y = "#00bfc4",
			high_xy = "#6a3d9a"
		),
		bivariateGreenRed = c(
			low_xy = "#e6e6e6",
			high_x = "#2ca25f",
			high_y = "#de2d26",
			high_xy = "#ffd92f"
		),
		bivariateBlueYellow = c(
			low_xy = "#e6e6e6",
			high_x = "#3182bd",
			high_y = "#fdd835",
			high_xy = "#31a354"
		)
	))
}

.bivariate_palette_names <- function() {
	return(names(.bivariate_palette_registry()))
}

.bivariate_palette_corners <- function(palette = "bivariateBright") {
	builtins <- .bivariate_palette_registry()
	corner_names <- c("low_xy", "high_x", "high_y", "high_xy")
	if (is.character(palette) && length(palette) == 1L && !is.na(palette)) {
		if (palette %in% names(builtins)) {
			return(builtins[[palette]])
		}
		stop(
			"`fg_palette` must be a built-in bivariate palette name or a vector of four colours when `fg_type = \"bivariate\"`.",
			call. = FALSE
		)
	}
	if (!is.character(palette) || length(palette) != 4L) {
		stop(
			"`fg_palette` must be a built-in bivariate palette name or a vector of four colours when `fg_type = \"bivariate\"`.",
			call. = FALSE
		)
	}
	if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
		if (!setequal(names(palette), corner_names)) {
			stop("Named custom bivariate palettes must contain colours named low_xy, high_x, high_y and high_xy.", call. = FALSE)
		}
		palette <- palette[corner_names]
	} else {
		names(palette) <- corner_names
	}
	if (!.is_colour_vector(unname(palette))) {
		stop(
			"`fg_palette` must be a built-in bivariate palette name or a vector of four colours when `fg_type = \"bivariate\"`.",
			call. = FALSE
		)
	}
	return(palette)
}

# Helper to create an n x n bivariate palette matrix by RGB bilinear blending.
.bivariate_palette <- function(palette = "bivariateBright", n = 9L) {
	if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2L) {
		stop("`n` must be a single numeric value greater than or equal to 2.", call. = FALSE)
	}
	n <- as.integer(n)
	corners <- grDevices::col2rgb(.bivariate_palette_corners(palette)) / 255
	out <- matrix(NA_character_, nrow = n, ncol = n)
	seq01 <- seq(0, 1, length.out = n)
	for (iy in seq_len(n)) {
		for (ix in seq_len(n)) {
			tx <- seq01[ix]
			ty <- seq01[iy]
			rgb <- (1 - tx) * (1 - ty) * corners[, "low_xy"] +
				tx * (1 - ty) * corners[, "high_x"] +
				(1 - tx) * ty * corners[, "high_y"] +
				tx * ty * corners[, "high_xy"]
			out[iy, ix] <- grDevices::rgb(rgb[1], rgb[2], rgb[3])
		}
	}
	return(out)
}

# Helper to resolve global categorical values for consistent panel palettes.
.resolve_categorical_values <- function(values, fallback) {
	if (is.null(values)) {
		values <- unique(na.omit(as.character(fallback)))
	} else {
		values <- unique(na.omit(as.character(values)))
	}
	return(values)
}

# Helper to compute shared limits and categorical values before panel splitting.
.resolve_grouped_layer_settings <- function(sce,
																					 background,
																					 foreground,
																					 bg_type,
																					 fg_type,
																					 bg_assay,
																					 fg_assay,
																					 fg_threshold,
																					 subset_idx,
																					 bg_limits,
																					 fg_limits) {
	bg <- .resolve_overlay_variable(
		sce = sce,
		variable = background,
		layer_name = "background",
		assay = bg_assay,
		type = bg_type
	)
	if (identical(fg_type, "bivariate")) {
		fg <- .prepare_bivariate_values(
			sce = sce,
			features = foreground,
			assay = fg_assay,
			limits = NULL
		)
	} else if (identical(fg_type, "combination")) {
		fg <- .prepare_combination_values(
			sce = sce,
			features = foreground,
			assay = fg_assay,
			threshold = fg_threshold
		)
	} else {
		fg <- .resolve_overlay_variable(
			sce = sce,
			variable = foreground,
			layer_name = "foreground",
			assay = fg_assay,
			type = fg_type
		)
	}

	resolved_bg_limits <- .resolve_layer_limits(
		limits = bg_limits,
		values = bg$values,
		type = bg$type,
		arg_name = "bg_limits",
		layer_name = "background"
	)
	fg_values_for_limits <- if (identical(fg$type, "none")) {
		fg$values
	} else if (identical(fg$type, "bivariate")) {
		NULL
	} else if (identical(fg$type, "combination")) {
		NULL
	} else {
		fg$values[subset_idx]
	}
	if (identical(fg$type, "bivariate")) {
		if (is.character(fg_limits) && identical(fg_limits, "shared")) {
			resolved_fg_limits <- .new_bivariate_limits(
				.range_or_null(fg$values_x[subset_idx], "fg_limits$x"),
				.range_or_null(fg$values_y[subset_idx], "fg_limits$y")
			)
		} else if (inherits(fg_limits, "scOverlay_bivariate_limits")) {
			resolved_fg_limits <- fg_limits
		} else if (is.null(fg_limits)) {
			resolved_fg_limits <- NULL
		} else {
			.validate_bivariate_fg_limits(fg_limits)
		}
	} else if (identical(fg$type, "combination")) {
		.validate_combination_fg_limits(fg_limits)
		resolved_fg_limits <- NULL
	} else {
		resolved_fg_limits <- .resolve_layer_limits(
			limits = fg_limits,
			values = fg_values_for_limits,
			type = fg$type,
			arg_name = "fg_limits",
			layer_name = "foreground"
		)
	}

	bg_values <- NULL
	if (identical(bg$type, "categorical")) {
		bg_values <- .resolve_categorical_values(NULL, bg$values)
	}
	fg_values <- NULL
	if (identical(fg$type, "categorical")) {
		fg_values <- .resolve_categorical_values(NULL, fg$values[subset_idx])
	} else if (identical(fg$type, "combination")) {
		fg_values <- levels(fg$labels)
	}

	return(list(
		bg_limits = resolved_bg_limits,
		fg_limits = resolved_fg_limits,
		bg_values = bg_values,
		fg_values = fg_values
	))
}

# Helper to validate rasterization controls
.validate_raster_args <- function(bg_raster, fg_raster, raster_dpi) {
	if (!.is_flag(bg_raster)) {
		stop("`bg_raster` must be TRUE or FALSE.", call. = FALSE)
	}
	if (!.is_flag(fg_raster)) {
		stop("`fg_raster` must be TRUE or FALSE.", call. = FALSE)
	}
	if (!is.numeric(raster_dpi) || length(raster_dpi) != 1L ||
			is.na(raster_dpi) || !is.finite(raster_dpi) || raster_dpi <= 0) {
		stop("`raster_dpi` must be a single positive finite numeric value.", call. = FALSE)
	}
	return(invisible(TRUE))
}

# Helper to construct vector or raster point layers
.geom_point_overlay <- function(data = NULL, mapping = NULL, raster = FALSE, raster_dpi = 300, ...) {
	if (isTRUE(raster)) {
		return(ggrastr::geom_point_rast(
			data = data,
			mapping = mapping,
			raster.dpi = raster_dpi,
			...
		))
	}

	return(ggplot2::geom_point(
		data = data,
		mapping = mapping,
		...
	))
}

# Helper to validate SingleCellExperiment inputs
.validate_sce <- function(sce) {
	if (!inherits(sce, "SingleCellExperiment")) {
		stop("`sce` must be a SingleCellExperiment object.", call. = FALSE)
	}
	return(invisible(TRUE))
}

# Helper to validate colData column arguments
.validate_coldata_column <- function(sce, column, arg_name) {
	if (!is.character(column) || length(column) != 1L || is.na(column)) {
		stop("`", arg_name, "` must be a single non-NA character string.", call. = FALSE)
	}

	cd_names <- colnames(SingleCellExperiment::colData(sce))
	if (!(column %in% cd_names)) {
		stop(
			"`", arg_name, "` must be one of colnames(colData(sce)): ",
			paste(cd_names, collapse = ", "),
			".",
			call. = FALSE
		)
	}

	return(invisible(TRUE))
}

# Helper to resolve a feature, colData column, or optional solid layer.
.resolve_overlay_variable <- function(sce,
																			variable,
																			layer_name,
																			assay = "logcounts",
																			type = c("auto", "continuous", "categorical")) {
	type <- match.arg(type)
	if (identical(variable, "solid")) {
		return(list(
			values = rep("all", ncol(sce)),
			type = "solid",
			source = "solid",
			label = variable
		))
	}
	if (identical(variable, "none")) {
		return(list(
			values = rep(NA_character_, ncol(sce)),
			type = "none",
			source = "none",
			label = variable
		))
	}

	if (variable %in% SummarizedExperiment::rownames(sce)) {
		assay_arg <- .layer_arg_name(layer_name, "assay")
		.validate_string_arg(assay, assay_arg)
		if (!(assay %in% SummarizedExperiment::assayNames(sce))) {
			stop(
				"`", assay_arg, "` must be one of assayNames(sce). Unknown assay: ",
				assay,
				".",
				call. = FALSE
			)
		}
		values <- as.numeric(SummarizedExperiment::assay(sce, assay)[variable, ])
		if (identical(type, "auto")) {
			type <- "continuous"
		}
		.validate_overlay_variable_type(values, type, layer_name)
		return(list(
			values = values,
			type = type,
			source = "assay",
			label = variable
		))
	}

	if (variable %in% colnames(SummarizedExperiment::colData(sce))) {
		values <- SummarizedExperiment::colData(sce)[[variable]]
		if (identical(type, "auto")) {
			type <- if (is.numeric(values)) "continuous" else "categorical"
		}
		.validate_overlay_variable_type(values, type, layer_name)
		return(list(
			values = values,
			type = type,
			source = "colData",
			label = variable
		))
	}

	stop(
		"`", layer_name, "` must be \"solid\", \"none\", a feature in rownames(sce), or a column in colData(sce). ",
		"Unknown ", layer_name, ": ", variable, ".",
		call. = FALSE
	)
}

# Helper to validate layer values against the requested scale type.
.validate_overlay_variable_type <- function(values, type, layer_name) {
	if (identical(type, "continuous") && !is.numeric(values)) {
		type_arg <- .layer_arg_name(layer_name, "type")
		stop(
			"`", type_arg, " = 'continuous'` requires a numeric ",
			layer_name,
			" variable.",
			call. = FALSE
		)
	}
	return(invisible(TRUE))
}

# Helper to map layer names to public argument names.
.layer_arg_name <- function(layer_name, suffix) {
	prefix <- switch(
		layer_name,
		background = "bg",
		foreground = "fg",
		layer_name
	)
	return(paste0(prefix, "_", suffix))
}

# Helper to extract and validate reduced dimension coordinates
.get_reduced_dim_coords <- function(sce, reduced_dim) {
	rd_names <- SingleCellExperiment::reducedDimNames(sce)
	if (!(reduced_dim %in% rd_names)) {
		stop(
			"`reduced_dim` must be one of reducedDimNames(sce): ",
			paste(rd_names, collapse = ", "),
			call. = FALSE
		)
	}

	rd <- SingleCellExperiment::reducedDim(sce, reduced_dim)
	if (ncol(rd) < 2L) {
		stop(
			"Reduced dimension `", reduced_dim, "` must contain at least two columns.",
			call. = FALSE
		)
	}

	coords <- as.data.frame(rd[, seq_len(2L), drop = FALSE])
	names(coords) <- c("X", "Y")
	return(coords)
}

# Helper to resolve the fixed colour used for background = "solid"
.resolve_solid_bg_colour <- function(bg_palette) {
	if (is.null(bg_palette)) return("gray70")
	if (!is.character(bg_palette) || length(bg_palette) != 1L || is.na(bg_palette)) {
		stop("`bg_palette` must be NULL or a single valid colour when `background = \"solid\"`.", call. = FALSE)
	}
	if (!.is_colour_vector(bg_palette)) {
		stop("`bg_palette` must be NULL or a single valid colour when `background = \"solid\"`.", call. = FALSE)
	}
	return(bg_palette)
}

# Helper to resolve the fixed colour used for foreground = "solid"
.resolve_solid_fg_colour <- function(fg_palette) {
	if (is.null(fg_palette)) return("red")
	if (!is.character(fg_palette) || length(fg_palette) != 1L || is.na(fg_palette)) {
		stop("`fg_palette` must be NULL or a single valid colour when `foreground = \"solid\"`.", call. = FALSE)
	}
	if (!.is_colour_vector(fg_palette)) {
		stop("`fg_palette` must be NULL or a single valid colour when `foreground = \"solid\"`.", call. = FALSE)
	}
	return(fg_palette)
}

# Helper to order foreground rows before plotting
.order_foreground <- function(coords_fg, fg_order) {
	fg_order <- match.arg(fg_order, c("input", "random", "ascending", "descending"))

	if (identical(fg_order, "input")) {
		return(coords_fg)
	}

	if (nrow(coords_fg) == 0L) {
		return(coords_fg)
	}

	if (identical(fg_order, "random")) {
		idx <- sample(seq_len(nrow(coords_fg)))
		return(coords_fg[idx, , drop = FALSE])
	}

	foreground <- coords_fg$foreground

	if (identical(fg_order, "ascending")) {
		idx <- order(foreground, na.last = FALSE)
		return(coords_fg[idx, , drop = FALSE])
	}

	if (identical(fg_order, "descending")) {
		idx <- order(foreground, decreasing = TRUE, na.last = FALSE)
		return(coords_fg[idx, , drop = FALSE])
	}

	stop("`fg_order` must be one of \"input\", \"random\", \"ascending\", or \"descending\".", call. = FALSE)
}

# Helper to avoid repeating the expected empty-subset warning in grouped plots.
.with_empty_subset_warning_muffled <- function(expr) {
	msg <- "`fg_subset_cells` did not match any cells; drawing an empty foreground layer."
	value <- withCallingHandlers(
		expr,
		warning = function(w) {
			if (identical(conditionMessage(w), msg)) {
				invokeRestart("muffleWarning")
			}
		}
	)
	return(value)
}

# Helper to evaluate fg_subset_cells
.eval_subset <- function(fg_subset_cells, sce) {
	ncells <- ncol(sce)

	validate_subset <- function(x, source) {
		if (!is.logical(x)) {
			stop("`fg_subset_cells` evaluated from ", source, " must return a logical vector.", call. = FALSE)
		}
		if (length(x) != ncells) {
			stop("`fg_subset_cells` evaluated from ", source, " must have length ncol(sce) (", ncells, "), not ", length(x), ".", call. = FALSE)
		}
		if (anyNA(x)) {
			warning("`fg_subset_cells` contains NA values; treating them as FALSE.", call. = FALSE)
			x[is.na(x)] <- FALSE
		}
		if (!any(x)) {
			warning("`fg_subset_cells` did not match any cells; drawing an empty foreground layer.", call. = FALSE)
		}
		return(x)
	}

	if (is.null(fg_subset_cells)) {
		return(rep(TRUE, ncells))
	} else if (is.logical(fg_subset_cells)) {
		return(validate_subset(fg_subset_cells, "a logical vector"))
	} else if (is.character(fg_subset_cells)) {
		cell_names <- colnames(sce)
		if (is.null(cell_names)) {
			stop("`fg_subset_cells` cannot select cell names because colnames(sce) is NULL.", call. = FALSE)
		}
		missing_names <- setdiff(fg_subset_cells, cell_names)
		if (length(missing_names) > 0L) {
			stop("`fg_subset_cells` contains cell names not present in colnames(sce): ", paste(missing_names, collapse = ", "),	".", call. = FALSE)
		}
		return(validate_subset(cell_names %in% fg_subset_cells, "cell names"))
	} else if (is.function(fg_subset_cells)) {
		res <- fg_subset_cells(sce)
		return(validate_subset(res, "a function"))
	} else if (is.language(fg_subset_cells)) { # quoted expression
		res <- tryCatch(
			with(as.data.frame(colData(sce)), eval(fg_subset_cells)),
			error = function(e) {
				stop("`fg_subset_cells` expression could not be evaluated against colData(sce): ", conditionMessage(e),	call. = FALSE	)
			}
		)
		return(validate_subset(res, "an expression"))
	} else {
		stop("`fg_subset_cells` must be NULL, a logical vector, a character vector of cell names, a function, or an expression.",	call. = FALSE)
	}
}
