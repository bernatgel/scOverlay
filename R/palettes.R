# Internal palette registry. Keep this object private; expose palette lookup
# through getPalette() and listPalettes().
scOverlay_palettes <- list(
	categorical = list(
		scOverlay = c(
			"1" = "#05d165",
			"2" = "#e6194b",
			"3" = "#ebad00",
			"4" = "#be6af3",
			"5" = "deepskyblue",
			"6" = "#fa70ff",
			"7" = "#46f0f0",
			"8" = "#ffe119",
			"9" = "#bcf60c",
			"10" = "#e6beff",
			"11" = "#008080",
			"12" = "#8baca3",
			"13" = "#9a6324",
			"14" = "#3688c9",
			"15" = "#800000",
			"16" = "#ff8270",
			"17" = "#808000",
			"18" = "#ffd8b1",
			"19" = "#000075",
			"20" = "#808080",
			"21" = "#fffac8",
			"22" = "orangered",
			"23" = "orange"
		),
		base = c(
			"1" = "#3B7EA1",
			"2" = "#D95F02",
			"3" = "#4DAF4A",
			"4" = "#984EA3",
			"5" = "#E41A1C",
			"6" = "#A6761D",
			"7" = "#1B9E77",
			"8" = "#E7298A",
			"9" = "#66A61E",
			"10" = "#7570B3",
			"11" = "#E6AB02",
			"12" = "#A6CEE3",
			"13" = "#B15928",
			"14" = "#F781BF",
			"15" = "#33A02C",
			"16" = "#6A3D9A",
			"17" = "#FF7F00",
			"18" = "#8DD3C7",
			"19" = "#FB9A99",
			"20" = "#BC80BD",
			"21" = "#FDB462",
			"22" = "#80B1D3",
			"23" = "#B3DE69",
			"24" = "#FCCDE5"
		),
		dark = c(
			"1" = "#1F4E79",
			"2" = "#B4463A",
			"3" = "#4C8C4A",
			"4" = "#7B5AA6",
			"5" = "#D98C21",
			"6" = "#2A9D8F",
			"7" = "#C7528B",
			"8" = "#6B6E23",
			"9" = "#8C564B",
			"10" = "#3A6EA5",
			"11" = "#E0B02E",
			"12" = "#5F9EA0",
			"13" = "#A23E48",
			"14" = "#6A994E",
			"15" = "#9A6FB0",
			"16" = "#C97B63",
			"17" = "#2F6F73",
			"18" = "#C77DFF",
			"19" = "#8D99AE",
			"20" = "#BC6C25",
			"21" = "#457B9D",
			"22" = "#E76F51",
			"23" = "#588157",
			"24" = "#9D4EDD"
		),
		soft = c(
			"1" = "#4E79A7",
			"2" = "#F28E2B",
			"3" = "#59A14F",
			"4" = "#B07AA1",
			"5" = "#E15759",
			"6" = "#76B7B2",
			"7" = "#EDC948",
			"8" = "#9C755F",
			"9" = "#FF9DA7",
			"10" = "#BAB0AC",
			"11" = "#86BCB6",
			"12" = "#D37295",
			"13" = "#8CD17D",
			"14" = "#B6992D",
			"15" = "#A0CBE8",
			"16" = "#FFBE7D",
			"17" = "#C5B0D5",
			"18" = "#FABFD2",
			"19" = "#79706E",
			"20" = "#D4A6C8",
			"21" = "#6B9AC4",
			"22" = "#E89C6C",
			"23" = "#88B04B",
			"24" = "#A26769"
		),
		bright = c(
			"1" = "#0057FF",
			"2" = "#FF1F1F",
			"3" = "#00A651",
			"4" = "#8A2BE2",
			"5" = "#FF8C00",
			"6" = "#00BFC4",
			"7" = "#FF00A8",
			"8" = "#7CAE00",
			"9" = "#FFD400",
			"10" = "#7F3C00",
			"11" = "#00E5FF",
			"12" = "#C000FF",
			"13" = "#FF5A00",
			"14" = "#008B8B",
			"15" = "#FF6EC7",
			"16" = "#4B0082",
			"17" = "#39FF14",
			"18" = "#B22222",
			"19" = "#1E90FF",
			"20" = "#9ACD32",
			"21" = "#FF1493",
			"22" = "#00FA9A",
			"23" = "#9400D3",
			"24" = "#FFB000"
		),
		pastel = c(
			"1" = "#8ECAE6",
			"2" = "#FFB703",
			"3" = "#A7C957",
			"4" = "#CDB4DB",
			"5" = "#FFADAD",
			"6" = "#BDE0FE",
			"7" = "#B8E0D2",
			"8" = "#FFD6A5",
			"9" = "#FFC8DD",
			"10" = "#D0F4DE",
			"11" = "#A9DEF9",
			"12" = "#E4C1F9",
			"13" = "#FCF6BD",
			"14" = "#CDEAC0",
			"15" = "#F1C0E8",
			"16" = "#CFBAF0",
			"17" = "#A3C4F3",
			"18" = "#90DBF4",
			"19" = "#98F5E1",
			"20" = "#FBF8CC",
			"21" = "#FDE4CF",
			"22" = "#FFCFD2",
			"23" = "#C1FBA4",
			"24" = "#B9FBC0"
		),
		combinationBright = c(
			"#D9D9D9",
			"#ff4346",
			"#55b2ff",
			"#68d865",
			"#984EA3",
			"#FF7F00",
			"#FFFF33",
			"#A65628",
			"#F781BF",
			"#00CED1",
			"#B3DE69",
			"#FDB462",
			"#80B1D3",
			"#FB8072",
			"#BC80BD",
			"#8DD3C7",
			"#CCEBC5",
			"#FFED6F",
			"#1B9E77",
			"#D95F02",
			"#7570B3",
			"#E7298A",
			"#66A61E",
			"#E6AB02",
			"#A6761D",
			"#666666",
			"#1F78B4",
			"#33A02C",
			"#FB9A99",
			"#CAB2D6",
			"#FDBF6F",
			"#B2DF8A"
		)
	),
	continuous = list(
		magma = c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF"),
		gray_blue = c("#c9c9c9", "#C7DCEF", "#7FB3D5", "#2874A6", "#08306B"),
		gray_red = c("#c9c9c9", "#F4B6B6", "#EF6F6C", "#CB181D", "#67000D"),
		gray_purple = c("#c9c9c9", "#D7BDE2", "#AF7AC5", "#7D3C98", "#3F007D"),
		blue_gold = c("#081D58", "#225EA8", "#41B6C4", "#C7E9B4", "#FFFFCC"),
		purple_orange = c("#2D004B", "#54278F", "#807DBA", "#FDB863", "#FFF7BC"),
		black_red_yellow = c("#000000", "#4A0000", "#B30000", "#F46D43", "#FFD92F", "#FFFFB2"),
		diverging_blue_white_red = c("#053061", "#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B", "#67001F"),
		diverging_green_white_purple = c("#00441B", "#238B45", "#A1D99B", "#F7F7F7", "#BCBDDC", "#756BB1", "#3F007D"),
		diverging_teal_gray_magenta = c("#004D40", "#00897B", "#80CBC4", "#F0F0F0", "#F4A3C1", "#D81B60", "#880E4F")
	)
)

#' Resolve a scOverlay palette
#'
#' \code{getPalette} converts a palette specification into a vector of colours
#' that can be used by scOverlay or directly in a \pkg{ggplot2} scale.
#'
#' Palette specifications can be \code{NULL}, the name of one of the palettes
#' bundled with scOverlay, a vector of colours, or a palette function. Palette
#' names are first looked up in scOverlay's internal palette registry for the
#' requested \code{type}, then in \pkg{viridisLite}, and then among ColorBrewer
#' palettes through \pkg{scales}.
#'
#' @details
#' When \code{palette = NULL}, \code{getPalette} uses the default scOverlay
#' palette for the requested type: \code{"scOverlay"} for categorical palettes and
#' \code{"gray_red"} for continuous palettes.
#'
#' Continuous palettes are represented as anchor colours and interpolated
#' internally to the requested number of colours. This means that a user supplied
#' vector such as \code{c("gray90", "red", "darkred")} can be used as a
#' continuous gradient.
#'
#' Categorical palettes can be named. When \code{values} is supplied and all
#' values are present in the palette names, colours are matched by name. When no
#' values are present in the palette names, colours are assigned by position. If
#' only some values match, \code{getPalette()} raises an error to avoid silent
#' colour remapping. Built-in categorical palettes intentionally keep
#' integer-like names such as \code{"1"}, \code{"2"} and \code{"3"} because
#' these are useful for cluster labels. Use \code{unname()} to force positional
#' matching for a named palette, or \code{setNames()} to assign or replace names.
#'
#' If a categorical palette has fewer colours than requested, additional colours
#' are generated by interpolation.
#' The categorical palette \code{"combinationBright"} contains 32 colours and is
#' used by default for \code{fg_type = "combination"} foregrounds in
#' \code{\link{plotOverlay}}.
#'
#' @param palette Palette specification. It can be \code{NULL}, a single palette
#' name, a vector of colours or a palette function. A palette function should
#' take the number of requested colours as its first argument and return a colour
#' vector. Defaults to \code{NULL}.
#'
#' @param type Character value indicating the type of palette to return.
#' Accepted values are \code{"continuous"} and \code{"categorical"}. Defaults to
#' \code{"continuous"}.
#'
#' @param n Number of colours to return. If \code{NULL}, it is inferred from
#' \code{values} when available, otherwise a default of \code{256} is used.
#' Non-integer values of \code{n} are truncated with \code{as.integer()}.
#' Defaults to \code{NULL}.
#'
#' @param values Optional character vector with the expected categorical values.
#' When supplied and the resolved palette is named, names are used to map colours
#' to values. Defaults to \code{NULL}.
#'
#' @return A character vector of colours. For categorical palettes, the returned
#' vector is named with \code{values} when \code{values} is supplied.
#'
#' @seealso \code{\link{listPalettes}}, \code{\link{plotOverlay}}.
#'
#' @examples
#' getPalette("gray_red", type = "continuous", n = 5)
#'
#' getPalette(
#'     palette = c("#E6E6E6", "#EF6F6C", "#67000D"),
#'     type = "continuous",
#'     n = 5
#' )
#'
#' getPalette(
#'     palette = "scOverlay",
#'     type = "categorical",
#'     n = 3,
#'     values = c("A", "B", "C")
#' )
#'
#' getPalette("combinationBright", type = "categorical", n = 8)
#'
#' getPalette(
#'     palette = c(A = "red", B = "blue"),
#'     type = "categorical",
#'     values = c("B", "A")
#' )
#'
#' pal <- getPalette("scOverlay", type = "categorical", n = 4)
#' setNames(unname(pal), c("Nerve", "PNF", "ANF", "MPNST"))
#'
#' @export

getPalette <- function(palette = NULL,
											 type = c("continuous", "categorical"),
											 n = NULL,
											 values = NULL) {
	type <- match.arg(type)
	if (!is.null(values)) values <- as.character(values)
	n <- .palette_n(n, values)

	if (is.null(palette)) {
		palette <- .default_palette_name(type)
	}

	if (is.function(palette)) {
		cols <- tryCatch(palette(n), error = function(e) {
			stop("`palette` function failed: ", conditionMessage(e), call. = FALSE)
		})
		return(.finalize_palette(cols, type = type, n = n, values = values))
	}

	if (is.character(palette)) {
		if (length(palette) == 1L && !.is_colour_vector(palette)) {
			if (is.na(palette)) {
				stop("Unknown palette `NA`.", call. = FALSE)
			}
			cols <- .palette_from_name(palette, type = type, n = n)
			return(.finalize_palette(cols, type = type, n = n, values = values))
		}
		return(.finalize_palette(palette, type = type, n = n, values = values))
	}

	stop("`palette` must be NULL, a character palette name, a colour vector, or a function.", call. = FALSE)
}



#' List available scOverlay palettes
#'
#' \code{listPalettes} returns metadata for the palettes bundled with scOverlay.
#'
#' @details
#' Categorical palettes return their stored colour vectors. Internally defined
#' continuous palettes return their stored anchor colours. Externally generated
#' continuous palettes, such as \pkg{viridisLite} palettes, return
#' \code{n = NA_integer_} and \code{colours = NA} because scOverlay does not
#' store a fixed colour vector for them. Bivariate palettes return their four
#' named corner colours: \code{low_xy}, \code{high_x}, \code{high_y} and
#' \code{high_xy}.
#'
#' The \code{is_default} column marks the palettes used by
#' \code{\link{getPalette}} when \code{palette = NULL}: \code{"scOverlay"} for
#' categorical palettes, \code{"gray_red"} for continuous palettes and
#' \code{"bivariateBright"} for bivariate foregrounds.
#'
#' @param type Character value indicating which palettes to list. Accepted
#' values are \code{"all"}, \code{"categorical"}, \code{"continuous"} and
#' \code{"bivariate"}. Defaults to \code{"all"}.
#'
#' @return A \code{data.frame} with palette \code{name}, \code{type}, stored
#' colour count \code{n}, a \code{colours} list-column containing stored palette
#' definitions, and a logical \code{is_default} column indicating the default
#' palettes.
#'
#' @seealso \code{\link{getPalette}}, \code{\link{plotPalettes}},
#' \code{\link{plotBivariateLegend}} and \code{\link{plotOverlay}}.
#'
#' @examples
#' listPalettes()
#'
#' listPalettes(type = "bivariate")
#'
#' pal_info <- listPalettes(type = "bivariate")
#' pal_info$colours[[which(pal_info$name == "bivariateBright")]]
#'
#' pal_info <- listPalettes(type = "continuous")
#' pal_info$colours[[which(pal_info$name == "gray_blue")]]
#'
#' @export

listPalettes <- function(type = c("all", "categorical", "continuous", "bivariate")) {
	type <- .palette_type(type)
	return(.palette_metadata(type))
}

#' Plot available scOverlay palettes
#'
#' \code{plotPalettes} displays the palettes bundled with scOverlay. Categorical
#' and continuous palettes are displayed as one-dimensional swatches. Bivariate
#' palettes are displayed as two-dimensional colour grids using the same
#' orientation as \code{\link{plotBivariateLegend}}: feature 1 increases from
#' left to right and feature 2 increases from bottom to top.
#'
#' Continuous palettes may be sampled or interpolated for display, including
#' externally generated palettes such as \pkg{viridisLite} palettes. Use
#' \code{\link{listPalettes}} to inspect, copy or modify stored palette
#' definitions.
#'
#' @param type Character value indicating which palettes to display. Accepted
#' values are \code{"all"}, \code{"categorical"}, \code{"continuous"} and
#' \code{"bivariate"}. Defaults to \code{"all"}.
#'
#' @param n Number of colours used to display continuous palettes, the maximum
#' number of colours shown for categorical palettes, and the number of bins per
#' axis for bivariate palettes. Non-integer values are truncated with
#' \code{as.integer()}. Defaults to \code{24}.
#'
#' @param show_names Logical value. If \code{TRUE}, categorical palette previews
#' show stored colour names when available, or position labels otherwise.
#' Bivariate palette names are always shown as plot titles. Defaults to
#' \code{TRUE}.
#'
#' @return A \pkg{ggplot2} object for a single preview section or a
#' \pkg{patchwork} object when multiple sections are shown.
#'
#' @seealso \code{\link{listPalettes}}, \code{\link{plotBivariateLegend}} and
#' \code{\link{getPalette}}.
#'
#' @examples
#' plotPalettes()
#' plotPalettes(type = "continuous")
#' plotPalettes(type = "bivariate")
#'
#' @export

plotPalettes <- function(type = c("all", "categorical", "continuous", "bivariate"),
												 n = 24,
												 show_names = TRUE) {
	type <- .palette_type(type)
	n <- .palette_n(n)
	if (!is.logical(show_names) || length(show_names) != 1L || is.na(show_names)) {
		stop("`show_names` must be TRUE or FALSE.", call. = FALSE)
	}
	info <- .palette_metadata(type)
	plots <- list()
	univariate_info <- info[info$type %in% c("categorical", "continuous"), , drop = FALSE]
	if (nrow(univariate_info) > 0L) {
		plots[["univariate"]] <- .plot_univariate_palettes(univariate_info, n, show_names)
	}
	bivariate_info <- info[info$type == "bivariate", , drop = FALSE]
	if (nrow(bivariate_info) > 0L) {
		plots[["bivariate"]] <- .plot_bivariate_palettes(bivariate_info, n)
	}
	if (length(plots) == 1L) {
		return(plots[[1L]])
	}
	return(patchwork::wrap_plots(plots, ncol = 1L, heights = c(2, 1)))
}

.palette_metadata <- function(type) {
	types <- if (identical(type, "all")) c("categorical", "continuous", "bivariate") else type
	info <- do.call(rbind, lapply(types, function(palette_type) {
		if (identical(palette_type, "bivariate")) {
			palette_names <- .bivariate_palette_names()
			n_colours <- rep(4L, length(palette_names))
			colours <- lapply(palette_names, .bivariate_palette_corners)
		} else {
			palettes <- scOverlay_palettes[[palette_type]]
			palette_names <- names(palettes)
			n_colours <- vapply(palettes, length, integer(1))
			colours <- unname(palettes)
			if (identical(palette_type, "continuous")) {
				viridis_names <- setdiff(.viridis_palette_names(), palette_names)
				palette_names <- c(palette_names, viridis_names)
				n_colours <- c(n_colours, rep(NA_integer_, length(viridis_names)))
				colours <- c(colours, as.list(rep(NA_character_, length(viridis_names))))
			}
		}
		info <- data.frame(
			name = palette_names,
			type = palette_type,
			n = n_colours,
			is_default = palette_names == .default_palette_name(palette_type),
			stringsAsFactors = FALSE
		)
		info <- .transform_palette_colours(info, colours)
		return(info)
	}))
	rownames(info) <- NULL
	return(info)
}

.palette_type <- function(type) {
	choices <- c("all", "categorical", "continuous", "bivariate")
	if (identical(type, choices)) {
		return("all")
	}
	if (!is.character(type) || length(type) != 1L || is.na(type) || !(type %in% choices)) {
		stop(
			"`type` must be one of \"all\", \"categorical\", \"continuous\" or \"bivariate\".",
			call. = FALSE
		)
	}
	return(type)
}

.transform_palette_colours <- function(info, colours) {
	info$colours <- I(colours)
	info <- info[, c("name", "type", "n", "colours", "is_default"), drop = FALSE]
	return(info)
}

.plot_univariate_palettes <- function(info, n, show_names) {
	swatch_data <- do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
		palette_name <- info$name[[i]]
		palette_type <- info$type[[i]]
		palette_label <- paste(palette_type, palette_name, sep = ": ")
		if (identical(palette_type, "continuous")) {
			cols <- getPalette(palette_name, type = "continuous", n = n)
			labels <- rep(NA_character_, length(cols))
		} else {
			cols <- scOverlay_palettes$categorical[[palette_name]]
			cols <- cols[seq_len(min(length(cols), n))]
			labels <- names(cols)
			if (is.null(labels) || all(is.na(labels) | !nzchar(labels))) {
				labels <- as.character(seq_along(cols))
			}
		}
		data.frame(
			palette = factor(palette_label, levels = rev(paste(info$type, info$name, sep = ": "))),
			type = palette_type,
			index = seq_along(cols),
			colour = unname(cols),
			label = labels,
			show_label = identical(palette_type, "categorical") && isTRUE(show_names),
			stringsAsFactors = FALSE
		)
	}))
	rownames(swatch_data) <- NULL
	label_data <- swatch_data[swatch_data$show_label, , drop = FALSE]

	p <- ggplot2::ggplot(
		swatch_data,
		ggplot2::aes(x = .data$index, y = .data$palette, fill = .data$colour)
	) +
		ggplot2::geom_tile(width = 1, height = 0.9) +
		ggplot2::geom_text(
			data = label_data,
			ggplot2::aes(x = .data$index, y = .data$palette),
			label = label_data$label,
			inherit.aes = FALSE,
			size = 2.4,
			colour = "black"
		) +
		ggplot2::scale_fill_identity() +
		ggplot2::facet_grid(ggplot2::vars(.data$type), scales = "free_y", space = "free_y") +
		ggplot2::labs(x = NULL, y = NULL) +
		ggplot2::theme_minimal(base_size = 11) +
		ggplot2::theme(
			axis.text.x = ggplot2::element_blank(),
			panel.grid = ggplot2::element_blank(),
			strip.text.y = ggplot2::element_text(angle = 0)
		)

	return(p)
}

.plot_bivariate_palettes <- function(info, n) {
	plots <- lapply(info$name, function(palette_name) {
		.plot_bivariate_palette_grid(
			palette = palette_name,
			features = c("Feature 1", "Feature 2"),
			n_bins = max(2L, n),
			title = palette_name
		)
	})
	names(plots) <- info$name
	return(patchwork::wrap_plots(plots))
}

# Resolve the requested palette length. If n is omitted for a categorical palette,
# use the number of expected values; otherwise fall back to a smooth gradient size.
.palette_n <- function(n, values = NULL) {
	if (is.null(n)) {
		if (!is.null(values)) return(max(1L, length(values)))
		return(256L)
	}
	n <- as.integer(n)
	if (length(n) != 1L || is.na(n) || n < 1L) {
		stop("`n` must be a single positive integer.", call. = FALSE)
	}
	return(n)
}

# Resolve the internal palette name used when users request package defaults.
.default_palette_name <- function(type) {
	type <- match.arg(type, c("categorical", "continuous", "bivariate"))
	if (identical(type, "categorical")) {
		return("scOverlay")
	}
	if (identical(type, "continuous")) {
		return("gray_red")
	}
	if (identical(type, "bivariate")) {
		return("bivariateBright")
	}
	stop("Invalid palette type.", call. = FALSE)
}

# Resolve a single palette name. scOverlay palettes take precedence for the
# requested type, followed by viridisLite functions and then ColorBrewer palettes.
.palette_from_name <- function(name, type, n) {
	palettes <- scOverlay_palettes[[type]]
	if (name %in% names(palettes)) {
		return(palettes[[name]])
	}

	viridis <- .viridis_palette(name, n)
	if (!is.null(viridis)) return(viridis)

	brewer <- .brewer_palette(name, type = type, n = n)
	if (!is.null(brewer)) return(brewer)

	stop("Unknown palette `", name, "`.", call. = FALSE)
}

# Look up viridisLite palettes with explicit namespace references so static
# package checks can see the imported namespace is used.
.viridis_palette <- function(name, n) {
	fun <- switch(
		name,
		viridis = viridisLite::viridis,
		magma = viridisLite::magma,
		inferno = viridisLite::inferno,
		plasma = viridisLite::plasma,
		cividis = viridisLite::cividis,
		rocket = viridisLite::rocket,
		mako = viridisLite::mako,
		turbo = viridisLite::turbo,
		NULL
	)
	if (is.null(fun)) {
		return(NULL)
	}
	return(fun(n))
}

.viridis_palette_names <- function() {
	return(c(
		"viridis",
		"magma",
		"inferno",
		"plasma",
		"cividis",
		"rocket",
		"mako",
		"turbo"
	))
}

# Look up a ColorBrewer palette by exact name. The type argument is advisory;
# the registered Brewer category is used so names resolve for either output type.
.brewer_palette <- function(name, type, n) {
	brewer_info <- RColorBrewer::brewer.pal.info
	if (!(name %in% rownames(brewer_info))) {
		return(NULL)
	}

	palette_info <- brewer_info[name, , drop = FALSE]
	brewer_category <- as.character(palette_info$category[[1L]])
	minimum_colours <- switch(
		brewer_category,
		seq = 3L,
		qual = 3L,
		div = 3L,
		NULL
	)
	if (is.null(minimum_colours)) {
		return(NULL)
	}

	maximum_colours <- as.integer(palette_info$maxcolors[[1L]])
	base_n <- min(max(n, minimum_colours), maximum_colours)
	cols <- RColorBrewer::brewer.pal(base_n, name)
	if (n < minimum_colours) {
		if (identical(brewer_category, "qual")) {
			cols <- cols[seq_len(n)]
		} else if (n == 1L) {
			cols <- cols[[2L]]
		} else {
			cols <- cols[c(1L, 3L)]
		}
	}
	return(cols)
}

# Validate, resize, and name a resolved colour vector. This is where categorical
# named palettes are matched to values and short palettes are interpolated.
.finalize_palette <- function(cols, type, n, values = NULL) {
	if (!is.character(cols) || length(cols) < 1L) {
		stop("`palette` must resolve to one or more colours.", call. = FALSE)
	}
	if (!.is_colour_vector(unname(cols))) {
		stop("`palette` contains invalid colours.", call. = FALSE)
	}

	if (type == "continuous") {
		pal <- .interpolate_palette(unname(cols), n)
		return(pal)
	}

	if (!is.null(values) && length(values) > 0L && !is.null(names(cols))) {
		matched <- match(values, names(cols))
		has_match <- !is.na(matched)
		if (all(has_match)) {
			pal <- unname(cols[matched])
			names(pal) <- values
			return(pal)
		}
		if (any(has_match)) {
			missing_values <- unique(values[!has_match])
			stop(
				"`palette` is a named categorical palette but does not define colours for: ",
				paste(missing_values, collapse = ", "),
				".",
				call. = FALSE
			)
		}
	}

	if (length(cols) < n) {
		warning(
			"Categorical palette has fewer colours than requested; interpolating to ",
			n,
			" colours.",
			call. = FALSE
		)
		pal <- .interpolate_palette(unname(cols), n)
	} else {
		pal <- unname(cols[seq_len(n)])
	}

	if (!is.null(values) && length(values) == length(pal)) names(pal) <- values
	return(pal)
}

# Expand a colour vector to the requested length. Single-colour palettes are
# repeated; multi-colour palettes become smooth gradients.
.interpolate_palette <- function(cols, n) {
	if (length(cols) == 1L) return(rep(cols, n))
	return(grDevices::colorRampPalette(cols)(n))
}

# Return TRUE when every string is understood by R as a colour. This accepts
# standard colour names and hex values.
.is_colour_vector <- function(cols) {
	if (!is.character(cols) || length(cols) < 1L) return(FALSE)
	return(all(vapply(cols, function(col) {
		!inherits(try(grDevices::col2rgb(col), silent = TRUE), "try-error")
	}, logical(1))))
}
