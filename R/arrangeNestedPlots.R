#' Arrange plots into a nested patchwork layout
#'
#' Internal helper used by \code{\link{plotOverlayPerNestedGroup}} to arrange a
#' named list of plots in rows defined by an outer grouping variable. Each row
#' contains inner groups ordered according to \code{group_order} when provided.
#' Compact layout skips missing outer/inner combinations and right-pads shorter
#' rows with blank spacers so panel widths stay aligned. Grid layout inserts
#' blank spacers for missing combinations in a full rectangular grid.
#'
#' @param sce A \linkS4class[SingleCellExperiment]{SingleCellExperiment} object.
#' @param plots A named list of \pkg{ggplot2} or \pkg{patchwork} plots. Names
#' should correspond to inner group values, or preferably to nested keys of the
#' form \code{outer__group}.
#' @param group_col A single character value with the name of the
#' \code{colData(sce)} column defining the inner groups.
#' @param outer_col A single character value with the name of the
#' \code{colData(sce)} column defining the outer groups.
#' @param group_order Optional character vector specifying the order of the
#' inner groups. Defaults to the order observed in \code{group_col}.
#' @param outer_order Optional character vector specifying the order of the
#' outer groups. Defaults to the order observed in \code{outer_col}.
#' @param spacer Ignored. Kept only for internal compatibility.
#' @param layout Character value controlling the nested arrangement. Accepted
#' values are \code{"compact"} and \code{"grid"}.
#'
#' @return A \pkg{patchwork} object with one row per outer group and a consistent
#' number of columns across rows.
#'
#' @keywords internal
  
arrangeNestedPlots <- function(
        sce,
        plots,
        group_col,
        outer_col,
        group_order = NULL,
        outer_order = NULL,
        spacer = NULL,
        layout = c("compact", "grid")
) {
    .validate_sce(sce)
    .validate_coldata_column(sce, group_col, "group_col")
    .validate_coldata_column(sce, outer_col, "outer_col")

    if (!is.list(plots)) {
        stop("`plots` must be a named list of plots.", call. = FALSE)
    }

    if (is.null(names(plots)) || any(is.na(names(plots))) ||
            any(names(plots) == "")) {
        stop("`plots` must be a named list.", call. = FALSE)
    }

    layout <- match.arg(layout)
    panel_layout <- .nested_panel_layout(
        sce = sce,
        group_col = group_col,
        outer_col = outer_col,
        group_order = group_order,
        outer_order = outer_order,
        layout = layout
    )
    if (is.null(spacer)) {
        spacer <- .empty_nested_panel()
    }

    get_plot <- function(outer, group, present) {
        nested_key <- paste(outer, group, sep = "__")

        if (nested_key %in% names(plots)) {
            return(plots[[nested_key]])
        }

        if (group %in% names(plots)) {
            return(plots[[group]])
        }
        if (identical(layout, "grid") && !isTRUE(present)) {
            return(spacer)
        }
        return(NULL)
    }

    row_plots <- lapply(unique(panel_layout$outer), function(outer_value) {
        row_layout <- panel_layout[panel_layout$outer == outer_value, , drop = FALSE]
        plots_for_row <- lapply(seq_len(nrow(row_layout)), function(i) {
            return(get_plot(
                outer = row_layout$outer[[i]],
                group = row_layout$group[[i]],
                present = row_layout$present[[i]]
            ))
        })
        plots_for_row <- Filter(Negate(is.null), plots_for_row)
        if (length(plots_for_row) == 0L) {
            return(NULL)
        }
        return(plots_for_row)
    })
    row_plots <- Filter(Negate(is.null), row_plots)

    if (length(row_plots) == 0L) {
        stop("No plots match the observed nested groups.", call. = FALSE)
    }

    ncol_layout <- max(vapply(row_plots, length, integer(1)))
    row_plots <- lapply(row_plots, function(plots_for_row) {
        n_missing <- ncol_layout - length(plots_for_row)
        if (n_missing > 0L) {
            plots_for_row <- c(plots_for_row, rep(list(spacer), n_missing))
        }
        return(plots_for_row)
    })
    flat_plots <- unlist(row_plots, recursive = FALSE)

    if (length(flat_plots) == 1L) {
        return(flat_plots[[1L]])
    }
    arranged <- patchwork::wrap_plots(flat_plots, ncol = ncol_layout)

    return(arranged)
}

.nested_panel_layout <- function(sce,
                                 group_col,
                                 outer_col,
                                 group_order = NULL,
                                 outer_order = NULL,
                                 layout = c("compact", "grid")) {
    .validate_sce(sce)
    .validate_coldata_column(sce, group_col, "group_col")
    .validate_coldata_column(sce, outer_col, "outer_col")
    layout <- match.arg(layout)

    group_vals <- SingleCellExperiment::colData(sce)[[group_col]]
    outer_vals <- SingleCellExperiment::colData(sce)[[outer_col]]
    group_vals_chr <- as.character(group_vals)
    outer_vals_chr <- as.character(outer_vals)

    valid <- !is.na(group_vals_chr) & nzchar(group_vals_chr) &
        !is.na(outer_vals_chr) & nzchar(outer_vals_chr)
    if (!any(valid)) {
        stop(
            "`group_col` and `outer_col` do not define any valid nested groups.",
            call. = FALSE
        )
    }

    if (!is.null(group_order)) {
        groups <- as.character(group_order)
    } else if (is.factor(group_vals)) {
        groups <- as.character(levels(group_vals))
    } else {
        groups <- unique(group_vals_chr[valid])
    }
    groups <- groups[!is.na(groups) & nzchar(groups)]
    groups <- groups[!duplicated(groups)]

    if (!is.null(outer_order)) {
        outers <- as.character(outer_order)
    } else if (is.factor(outer_vals)) {
        outers <- as.character(levels(outer_vals))
    } else {
        outers <- unique(outer_vals_chr[valid])
    }
    outers <- outers[!is.na(outers) & nzchar(outers)]
    outers <- outers[!duplicated(outers)]

    panel_indices <- .build_nested_panel_indices(
        outer_vals = outer_vals_chr,
        group_vals = group_vals_chr,
        outer_order = outers,
        group_order = groups,
        layout = layout
    )
    keys <- names(panel_indices)
    key_parts <- strsplit(keys, "__", fixed = TRUE)
    out <- data.frame(
        outer = vapply(key_parts, `[[`, character(1), 1L),
        group = vapply(key_parts, `[[`, character(1), 2L),
        key = keys,
        present = vapply(panel_indices, any, logical(1)),
        stringsAsFactors = FALSE
    )

    if (is.null(out) || nrow(out) == 0L) {
        stop("No observed outer/group combinations are available.", call. = FALSE)
    }
    rownames(out) <- NULL
    return(out)
}

.nested_panel_index <- function(outer_vals, group_vals, outer, group) {
    idx <- outer_vals == outer & group_vals == group
    idx[is.na(idx)] <- FALSE
    return(idx)
}

.resolve_nested_layout <- function(layout, outer_vals, group_vals) {
    layout <- match.arg(layout, c("auto", "compact", "grid"))
    if (!identical(layout, "auto")) {
        return(layout)
    }
    outer_vals <- as.character(outer_vals)
    group_vals <- as.character(group_vals)
    ok <- !is.na(outer_vals) & nzchar(outer_vals) &
        !is.na(group_vals) & nzchar(group_vals)
    if (!any(ok)) {
        return("compact")
    }
    tab <- table(outer_vals[ok], group_vals[ok])
    n_outer_per_group <- colSums(tab > 0)
    if (any(n_outer_per_group > 1L)) {
        return("grid")
    }
    return("compact")
}

.resolve_nested_panel_titles <- function(panel_titles) {
    panel_titles <- match.arg(
        panel_titles,
        c("auto", "outer_inner", "inner", "none")
    )

    if (identical(panel_titles, "auto")) {
        panel_titles <- "outer_inner"
    }

    return(panel_titles)
}

.nested_panel_title <- function(outer, group, panel_titles) {
    if (identical(panel_titles, "none")) {
        return(NULL)
    }

    if (identical(panel_titles, "inner")) {
        return(group)
    }

    if (identical(panel_titles, "outer_inner")) {
        title <- paste(group, outer, sep = " - ")
        return(title)
    }

    return(NULL)
}

.build_nested_panel_indices <- function(outer_vals,
                                        group_vals,
                                        outer_order,
                                        group_order,
                                        layout = c("compact", "grid")) {
    layout <- match.arg(layout)
    outer_vals <- as.character(outer_vals)
    group_vals <- as.character(group_vals)
    panels <- list()
    for (outer in outer_order) {
        for (group in group_order) {
            key <- paste(outer, group, sep = "__")
            idx <- .nested_panel_index(outer_vals, group_vals, outer, group)
            if (!any(idx) && identical(layout, "compact")) {
                next
            }
            panels[[key]] <- idx
        }
    }
    return(panels)
}

.empty_nested_panel <- function() {
    p <- ggplot2::ggplot() + ggplot2::theme_void()
    return(p)
}
