patchwork_plot_count <- function(p) {
  return(length(p$patches$plots) + 1L)
}

patchwork_child_plots <- function(p) {
  return(c(list(p), p$patches$plots))
}

grouped_foreground_layer_data <- function(p) {
  layer_data <- lapply(p$layers, function(layer) layer$data)
  foreground_layers <- Filter(
    function(x) is.data.frame(x) && "foreground" %in% colnames(x),
    layer_data
  )
  expect_true(length(foreground_layers) >= 1L)
  return(foreground_layers[[length(foreground_layers)]])
}

grouped_plot_data <- function(p) {
  expect_true(is.data.frame(p$data))
  expect_true(all(c("X", "Y", "background", "foreground") %in% colnames(p$data)))
  return(p$data)
}

expected_nested_panel_count <- function(sce, group_col, outer_col) {
  layout <- scOverlay:::.nested_panel_layout(
    sce = sce,
    group_col = group_col,
    outer_col = outer_col
  )
  return(nrow(layout))
}

expected_nested_ncol <- function(sce, group_col, outer_col) {
  layout <- scOverlay:::.nested_panel_layout(
    sce = sce,
    group_col = group_col,
    outer_col = outer_col
  )
  return(max(table(layout$outer)))
}

make_repeated_nested_sce <- function(missing_o2_g1 = FALSE) {
  outer <- c("O1", "O1", "O2", "O2")
  group <- c("G1", "G2", "G1", "G2")
  if (isTRUE(missing_o2_g1)) {
    outer <- outer[-3]
    group <- group[-3]
  }
  ncells <- length(outer)
  cells <- paste0("cell", seq_len(ncells))
  assays <- list(
    logcounts = matrix(
      seq_len(ncells),
      nrow = 1L,
      dimnames = list("Gene1", cells)
    )
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = data.frame(
      outer = outer,
      group = group,
      row.names = cells
    )
  )
  SingleCellExperiment::reducedDim(sce, "TSNE") <- matrix(
    seq_len(ncells * 2L),
    ncol = 2L,
    dimnames = list(cells, c("TSNE1", "TSNE2"))
  )
  return(sce)
}

test_that("plotOverlayPerGroup returns patchwork object", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerGroup(
      sce = sce,
      foreground = "SOX10"
    ),
    "group_col"
  )

  g <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "solid"
  )

  expect_s3_class(g, "patchwork")

  expect_error(
    plotOverlayPerGroup(
      sce = sce,
      group_col = "sample_type",
      foreground = "SOX10",
      out_file = tempfile(fileext = ".pdf")
    ),
    "unused argument"
  )
})

test_that("plotOverlayPerGroup returns NULL for all-NA groups", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$sample_type <- NA_character_

  expect_warning(
    g <- plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample_type",
      background = "solid"
    ),
    "No valid groups"
  )
  expect_null(g)
})

test_that("plotOverlayPerGroup validates group_col", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "missing_group",
      background = "solid"
    ),
    "`group_col` must be one of colnames\\(colData\\(sce\\)\\)"
  )

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = c("sample", "sample_type"),
      background = "solid"
    ),
    "`group_col` must be a single non-NA character string"
  )

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = NA_character_,
      background = "solid"
    ),
    "`group_col` must be a single non-NA character string"
  )
})

test_that("plotOverlayPerGroup forwards rasterization arguments", {
  sce <- get_example_sce()

  g <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "solid",
    bg_raster = TRUE,
    fg_raster = TRUE,
    raster_dpi = 600
  )

  expect_s3_class(g, "patchwork")
})

test_that("plotOverlayPerGroup forwards background type and assay", {
  sce <- get_example_sce()

  g_continuous <- plotOverlayPerGroup(
    sce,
    foreground = "cluster",
    group_col = "sample_type",
    background = "SOX10",
    bg_type = "continuous",
    bg_assay = "counts",
    fg_type = "categorical"
  )
  expect_s3_class(g_continuous, "patchwork")

  g_categorical <- plotOverlayPerGroup(
    sce,
    foreground = "S100B",
    group_col = "sample_type",
    background = "SOX10",
    bg_type = "categorical",
    bg_assay = "counts",
    bg_palette = unname(scOverlay:::scOverlay_palettes$categorical$scOverlay)
  )
  expect_s3_class(g_categorical, "patchwork")
})

test_that("plotOverlayPerGroup works with solid foreground", {
  sce <- get_example_sce()

  g <- plotOverlayPerGroup(
    sce,
    foreground = "solid",
    group_col = "sample_type",
    background = "cluster",
    fg_subset_cells = quote(sample_type != "Nerve"),
    fg_palette = "red"
  )

  expect_s3_class(g, "patchwork")
})

test_that("plotOverlayPerGroup works with none layers", {
  sce <- get_example_sce()

  g_bg_none <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "none"
  )
  expect_s3_class(g_bg_none, "patchwork")

  g_fg_none <- plotOverlayPerGroup(
    sce,
    foreground = "none",
    group_col = "sample_type",
    background = "sample_type"
  )
  expect_s3_class(g_fg_none, "patchwork")
})

test_that("plotOverlayPerGroup supports shared limits and legends", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_fg <- seq_len(ncol(sce))
  SummarizedExperiment::colData(sce)$numeric_bg <- rev(seq_len(ncol(sce)))
  subset_idx <- SummarizedExperiment::colData(sce)$sample_type != "Nerve"

  g_fg <- plotOverlayPerGroup(
    sce,
    foreground = "numeric_fg",
    group_col = "sample_type",
    background = "solid",
    fg_limits = "shared",
    fg_subset_cells = subset_idx
  )
  expect_s3_class(g_fg, "patchwork")
  fg_limits <- range(SummarizedExperiment::colData(sce)$numeric_fg[subset_idx])
  lapply(patchwork_child_plots(g_fg), function(p) {
    expect_equal(p$scales$scales[[1]]$limits, fg_limits)
  })

  SummarizedExperiment::colData(sce)$numeric_y <- rev(seq_len(ncol(sce)))
  expect_silent(
    g_bivar_default_legend <- plotOverlayPerGroup(
      sce,
      foreground = c("numeric_fg", "numeric_y"),
      group_col = "sample_type",
      background = "solid",
      fg_type = "bivariate",
      fg_limits = "shared",
      fg_subset_cells = subset_idx
    )
  )
  expect_s3_class(g_bivar_default_legend, "patchwork")

  expect_silent(
    g_bivar <- plotOverlayPerGroup(
      sce,
      foreground = c("numeric_fg", "numeric_y"),
      group_col = "sample_type",
      background = "solid",
      fg_type = "bivariate",
      fg_limits = "shared",
      fg_subset_cells = subset_idx,
      fg_legend = FALSE
    )
  )
  expect_s3_class(g_bivar, "patchwork")
  bivar_layers <- lapply(patchwork_child_plots(g_bivar), grouped_foreground_layer_data)
  bivar_layers <- Filter(function(x) nrow(x) > 0L, bivar_layers)
  expect_true(all(vapply(bivar_layers, function(x) ".fg_colour" %in% colnames(x), logical(1))))

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = c("numeric_fg", "numeric_y"),
      group_col = "sample_type",
      background = "solid",
      fg_type = "bivariate",
      fg_limits = c(0, 5),
      fg_legend = FALSE
    ),
    "Numeric foreground limits are not supported"
  )

  expect_warning(
    plotOverlayPerGroup(
      sce,
      foreground = c("numeric_fg", "numeric_y"),
      group_col = "sample_type",
      background = "solid",
      fg_type = "bivariate",
      fg_legend = TRUE
    ),
    "plotBivariateLegend"
  )

  g_bg <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "numeric_bg",
    bg_limits = "shared",
    shared_legend = TRUE
  )
  expect_s3_class(g_bg, "patchwork")
  bg_limits <- range(SummarizedExperiment::colData(sce)$numeric_bg)
  lapply(patchwork_child_plots(g_bg), function(p) {
    expect_equal(p$scales$scales[[1]]$limits, bg_limits)
  })

  g_unshared <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "numeric_bg",
    shared_legend = FALSE
  )
  expect_s3_class(g_unshared, "patchwork")
})

test_that("plotOverlayPerGroup uses global categorical values across panels", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$panel_cat <- rep(c("A", "B"), length.out = ncol(sce))
  pal <- c(A = "red", B = "blue", C = "green")

  g <- plotOverlayPerGroup(
    sce,
    foreground = "panel_cat",
    group_col = "sample_type",
    background = "solid",
    fg_type = "categorical",
    fg_palette = pal,
    fg_subset_cells = rep(TRUE, ncol(sce))
  )
  expect_s3_class(g, "patchwork")
  lapply(patchwork_child_plots(g), function(p) {
    scale <- p$scales$scales[[1]]
    expect_equal(scale$limits, c("A", "B"))
    expect_equal(scale$map(c("A", "B")), c("red", "blue"))
  })

  g_bg <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    background = "panel_cat",
    bg_type = "categorical",
    bg_palette = pal
  )
  expect_s3_class(g_bg, "patchwork")
  lapply(patchwork_child_plots(g_bg), function(p) {
    scale <- p$scales$scales[[1]]
    expect_equal(scale$limits, c("A", "B"))
    expect_equal(scale$map(c("A", "B")), c("red", "blue"))
  })

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "panel_cat",
      group_col = "sample_type",
      background = "solid",
      fg_type = "categorical",
      fg_palette = c(A = "red")
    ),
    "does not define colours for: B"
  )
})

test_that("plotOverlayPerGroup preserves panels for empty foreground subsets", {
  sce <- get_example_sce()
  group_order <- c("Nerve", "PNF", "ANNUBP", "MPNST")

  expect_warning(
    g_empty <- plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample_type",
      group_order = group_order,
      background = "sample_type",
      fg_subset_cells = quote(sample_type == "NOT_PRESENT")
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(g_empty, "patchwork")
  expect_equal(patchwork_plot_count(g_empty), length(group_order))
  expect_equal(g_empty$patches$layout$ncol, length(group_order))

  g_partial <- plotOverlayPerGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample_type",
    group_order = group_order,
    background = "sample_type",
    fg_subset_cells = quote(sample_type == "MPNST")
  )
  expect_s3_class(g_partial, "patchwork")
  expect_equal(patchwork_plot_count(g_partial), length(group_order))
  expect_equal(g_partial$patches$layout$ncol, length(group_order))
})

test_that("plotOverlayPerGroup validates rasterization arguments early", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample_type",
      background = "solid",
      bg_raster = NA
    ),
    "`bg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample_type",
      background = "solid",
      fg_raster = c(TRUE, FALSE)
    ),
    "`fg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample_type",
      background = "solid",
      raster_dpi = Inf
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
})

test_that("plotOverlayPerNestedGroup returns patchwork object", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerNestedGroup(
      sce = sce,
      foreground = "SOX10"
    ),
    "group_col"
  )

  expect_error(
    plotOverlayPerNestedGroup(
      sce = sce,
      group_col = "sample",
      foreground = "SOX10"
    ),
    "outer_col"
  )

  g <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    background = "solid"
  )

  expect_s3_class(g, "patchwork")

  expect_error(
    plotOverlayPerNestedGroup(
      sce = sce,
      group_col = "sample",
      outer_col = "sample_type",
      foreground = "SOX10",
      width = 3
    ),
    "unused argument"
  )
})

test_that("plotOverlayPerNestedGroup returns NULL for all-NA outer groups", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$sample_type <- NA_character_

  expect_warning(
    g <- plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "sample_type",
      background = "solid"
    ),
    "No valid outer groups"
  )
  expect_null(g)
})

test_that("plotOverlayPerNestedGroup validates grouping columns", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "missing_group",
      outer_col = "sample_type",
      background = "solid"
    ),
    "`group_col` must be one of colnames\\(colData\\(sce\\)\\)"
  )

  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "missing_outer",
      background = "solid"
    ),
    "`outer_col` must be one of colnames\\(colData\\(sce\\)\\)"
  )

  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = c("sample", "sample_type"),
      outer_col = "sample_type",
      background = "solid"
    ),
    "`group_col` must be a single non-NA character string"
  )

  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = NA_character_,
      background = "solid"
    ),
    "`outer_col` must be a single non-NA character string"
  )
})

test_that("plotOverlayPerNestedGroup forwards rasterization arguments", {
  sce <- get_example_sce()

  g <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    background = "solid",
    bg_raster = TRUE,
    fg_raster = TRUE,
    raster_dpi = 600
  )

  expect_s3_class(g, "patchwork")
})

test_that("plotOverlayPerNestedGroup forwards background type and assay", {
  sce <- get_example_sce()

  g_continuous <- plotOverlayPerNestedGroup(
    sce,
    foreground = "cluster",
    group_col = "sample",
    outer_col = "sample_type",
    background = "SOX10",
    bg_type = "continuous",
    bg_assay = "counts",
    fg_type = "categorical"
  )
  expect_s3_class(g_continuous, "patchwork")

  g_categorical <- plotOverlayPerNestedGroup(
    sce,
    foreground = "S100B",
    group_col = "sample",
    outer_col = "sample_type",
    background = "SOX10",
    bg_type = "categorical",
    bg_assay = "counts",
    bg_palette = unname(scOverlay:::scOverlay_palettes$categorical$scOverlay)
  )
  expect_s3_class(g_categorical, "patchwork")
})

test_that("plotOverlayPerNestedGroup works with solid foreground", {
  sce <- get_example_sce()

  g <- plotOverlayPerNestedGroup(
    sce,
    foreground = "solid",
    group_col = "sample",
    outer_col = "sample_type",
    background = "cluster",
    fg_subset_cells = quote(sample_type != "Nerve"),
    fg_palette = "#E41A1C"
  )

  expect_s3_class(g, "patchwork")
})

test_that("plotOverlayPerNestedGroup works with none layers", {
  sce <- get_example_sce()

  g_bg_none <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    background = "none"
  )
  expect_s3_class(g_bg_none, "patchwork")

  g_fg_none <- plotOverlayPerNestedGroup(
    sce,
    foreground = "none",
    group_col = "sample",
    outer_col = "sample_type",
    background = "sample_type"
  )
  expect_s3_class(g_fg_none, "patchwork")
})

test_that("plotOverlayPerNestedGroup supports shared limits and legends", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_fg <- seq_len(ncol(sce))
  SummarizedExperiment::colData(sce)$numeric_bg <- rev(seq_len(ncol(sce)))

  g_fg <- plotOverlayPerNestedGroup(
    sce,
    foreground = "numeric_fg",
    group_col = "sample",
    outer_col = "sample_type",
    background = "solid",
    fg_limits = "shared",
    shared_legend = TRUE
  )
  expect_s3_class(g_fg, "patchwork")
  fg_limits <- range(SummarizedExperiment::colData(sce)$numeric_fg)
  lapply(patchwork_child_plots(g_fg), function(p) {
    if (length(p$scales$scales) > 0L) {
      expect_equal(p$scales$scales[[1]]$limits, fg_limits)
    }
  })

  g_bg <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    background = "numeric_bg",
    bg_limits = "shared"
  )
  expect_s3_class(g_bg, "patchwork")
  bg_limits <- range(SummarizedExperiment::colData(sce)$numeric_bg)
  lapply(patchwork_child_plots(g_bg), function(p) {
    if (length(p$scales$scales) > 0L) {
      expect_equal(p$scales$scales[[1]]$limits, bg_limits)
    }
  })
})

test_that("grouped plots support combination foregrounds", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$combo_a <- seq_len(ncol(sce)) %% 3
  SummarizedExperiment::colData(sce)$combo_b <- rev(seq_len(ncol(sce))) %% 4

  g <- plotOverlayPerGroup(
    sce,
    foreground = c("combo_a", "combo_b"),
    group_col = "sample_type",
    background = "solid",
    fg_type = "combination",
    fg_threshold = 0,
    fg_limits = "shared",
    shared_legend = TRUE
  )
  expect_s3_class(g, "patchwork")

  expect_error(
    plotOverlayPerGroup(
      sce,
      foreground = c("combo_a", "combo_b"),
      group_col = "sample_type",
      background = "solid",
      fg_type = "combination",
      fg_limits = c(0, 5)
    ),
    "`fg_limits` is not supported"
  )
})

test_that("nested grouped plots support combination foregrounds", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$combo_a <- seq_len(ncol(sce)) %% 3
  SummarizedExperiment::colData(sce)$combo_b <- rev(seq_len(ncol(sce))) %% 4

  g <- plotOverlayPerNestedGroup(
    sce,
    foreground = c("combo_a", "combo_b"),
    group_col = "sample",
    outer_col = "sample_type",
    background = "solid",
    fg_type = "combination",
    fg_threshold = 0,
    fg_limits = "shared",
    shared_legend = TRUE
  )
  expect_s3_class(g, "patchwork")
})

test_that("nested panel layout uses observed outer/group combinations", {
  sce <- make_repeated_nested_sce()

  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "auto",
      outer_vals = c("O1", "O1", "O2", "O2"),
      group_vals = c("G1", "G2", "G3", "G4")
    ),
    "compact"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "auto",
      outer_vals = c("O1", "O1", "O2", "O2"),
      group_vals = c("G1", "G2", "G1", "G2")
    ),
    "grid"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "compact",
      outer_vals = c("O1", NA, "O2"),
      group_vals = c("G1", "G1", NA)
    ),
    "compact"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "grid",
      outer_vals = c("O1", NA, "O2"),
      group_vals = c("G1", "G1", NA)
    ),
    "grid"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "auto",
      outer_vals = c("O1", NA, "O2"),
      group_vals = c("G1", "G1", NA)
    ),
    "compact"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "auto",
      outer_vals = c("O1", "O1", "", NA, "O2"),
      group_vals = c("G1", "G2", "G1", "G2", "G3")
    ),
    "compact"
  )
  expect_identical(
    scOverlay:::.resolve_nested_layout(
      "auto",
      outer_vals = c("O1", "O2"),
      group_vals = c("G1", "G1")
    ),
    "grid"
  )

  layout <- scOverlay:::.nested_panel_layout(
    sce = sce,
    group_col = "group",
    outer_col = "outer",
    group_order = c("G1", "G2"),
    outer_order = c("O1", "O2"),
    layout = "compact"
  )

  expect_equal(layout$key, c("O1__G1", "O1__G2", "O2__G1", "O2__G2"))
  expect_equal(layout$outer, c("O1", "O1", "O2", "O2"))
  expect_equal(layout$group, c("G1", "G2", "G1", "G2"))

  outer_vals <- as.character(SummarizedExperiment::colData(sce)$outer)
  group_vals <- as.character(SummarizedExperiment::colData(sce)$group)
  expect_equal(
    which(scOverlay:::.nested_panel_index(outer_vals, group_vals, "O1", "G1")),
    1L
  )
  expect_equal(
    which(scOverlay:::.nested_panel_index(outer_vals, group_vals, "O2", "G1")),
    3L
  )

  panel_indices <- scOverlay:::.build_nested_panel_indices(
    outer_vals = outer_vals,
    group_vals = group_vals,
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "grid"
  )
  expect_equal(names(panel_indices), c("O1__G1", "O1__G2", "O2__G1", "O2__G2"))
  expect_equal(which(panel_indices[["O1__G1"]]), 1L)
  expect_equal(which(panel_indices[["O1__G2"]]), 2L)
  expect_equal(which(panel_indices[["O2__G1"]]), 3L)
  expect_equal(which(panel_indices[["O2__G2"]]), 4L)

  expect_identical(scOverlay:::.resolve_nested_panel_titles("auto"), "outer_inner")
  expect_identical(scOverlay:::.nested_panel_title("O1", "G1", "outer_inner"), "G1 - O1")
  expect_identical(scOverlay:::.nested_panel_title("O1", "G1", "inner"), "G1")
  expect_null(scOverlay:::.nested_panel_title("O1", "G1", "none"))
})

test_that("nested grouped plots handle repeated inner groups by outer/group pair", {
  sce <- make_repeated_nested_sce()

  g <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "grid"
  )

  expect_s3_class(g, "patchwork")
  expect_equal(
    expected_nested_panel_count(sce, "group", "outer"),
    4L
  )
})

test_that("nested grouped plots skip missing combinations in compact layout", {
  sce <- make_repeated_nested_sce(missing_o2_g1 = TRUE)

  layout <- scOverlay:::.nested_panel_layout(
    sce = sce,
    group_col = "group",
    outer_col = "outer",
    group_order = c("G1", "G2"),
    outer_order = c("O1", "O2"),
    layout = "compact"
  )
  expect_equal(layout$key, c("O1__G1", "O1__G2", "O2__G2"))
  expect_false("O2__G1" %in% layout$key)

  grid_layout <- scOverlay:::.nested_panel_layout(
    sce = sce,
    group_col = "group",
    outer_col = "outer",
    group_order = c("G1", "G2"),
    outer_order = c("O1", "O2"),
    layout = "grid"
  )
  expect_equal(grid_layout$key, c("O1__G1", "O1__G2", "O2__G1", "O2__G2"))
  expect_false(grid_layout$present[grid_layout$key == "O2__G1"])
  expect_true(grid_layout$present[grid_layout$key == "O2__G2"])

  g <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "compact"
  )
  expect_s3_class(g, "patchwork")
  expect_equal(patchwork_plot_count(g), 4L)
  expect_equal(g$patches$layout$ncol, 2L)

  g_grid <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "grid"
  )
  expect_s3_class(g_grid, "patchwork")
  expect_equal(patchwork_plot_count(g_grid), 4L)
  expect_equal(g_grid$patches$layout$ncol, 2L)
  grid_children <- patchwork_child_plots(g_grid)
  real_children <- Filter(function(p) length(p$layers) > 0L, grid_children)
  spacer_children <- Filter(function(p) length(p$layers) == 0L, grid_children)
  expect_equal(length(real_children), 3L)
  expect_equal(length(spacer_children), 1L)
  lapply(real_children, function(p) {
    expect_equal(nrow(grouped_plot_data(p)), ncol(sce))
    expect_equal(nrow(grouped_foreground_layer_data(p)), 1L)
  })
})

test_that("nested grouped plots keep full background and panel foreground subset", {
  sce <- make_repeated_nested_sce()

  g <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "grid"
  )

  children <- patchwork_child_plots(g)
  expect_equal(length(children), 4L)
  lapply(children, function(p) {
    expect_equal(nrow(grouped_plot_data(p)), ncol(sce))
    expect_equal(nrow(grouped_foreground_layer_data(p)), 1L)
  })

  fg_x <- vapply(children, function(p) {
    grouped_foreground_layer_data(p)$X[[1L]]
  }, numeric(1))
  expect_setequal(fg_x, c(1, 2, 3, 4))
})

test_that("nested grouped plots support simple panel titles", {
  sce <- make_repeated_nested_sce(missing_o2_g1 = TRUE)

  g_auto <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "grid"
  )
  children <- patchwork_child_plots(g_auto)
  real_children <- Filter(function(p) length(p$layers) > 0L, children)
  spacer_children <- Filter(function(p) length(p$layers) == 0L, children)
  expect_setequal(
    vapply(real_children, function(p) p$labels$title, character(1)),
    c("G1 - O1", "G2 - O1", "G2 - O2")
  )
  expect_equal(length(spacer_children), 1L)
  expect_null(spacer_children[[1L]]$labels$title)

  g_inner <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "compact",
    panel_titles = "inner"
  )
  real_inner <- Filter(function(p) length(p$layers) > 0L, patchwork_child_plots(g_inner))
  expect_setequal(
    vapply(real_inner, function(p) p$labels$title, character(1)),
    c("G1", "G2")
  )

  g_none <- plotOverlayPerNestedGroup(
    sce,
    outer_col = "outer",
    group_col = "group",
    foreground = "Gene1",
    background = "solid",
    outer_order = c("O1", "O2"),
    group_order = c("G1", "G2"),
    layout = "compact",
    panel_titles = "none"
  )
  real_none <- Filter(function(p) length(p$layers) > 0L, patchwork_child_plots(g_none))
  lapply(real_none, function(p) {
    expect_null(p$labels$title)
  })
})

test_that("plotOverlayPerNestedGroup handles shared limits with empty foreground panels", {
  sce <- get_example_sce()

  g <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    background = "solid",
    fg_subset_cells = quote(sample_type == "MPNST"),
    fg_limits = "shared"
  )
  expect_s3_class(g, "patchwork")
  expect_equal(expected_nested_panel_count(sce, "sample", "sample_type"), 6L)
})

test_that("plotOverlayPerNestedGroup preserves layout for empty foreground subsets", {
  sce <- get_example_sce()
  sample_order <- unique(as.character(SummarizedExperiment::colData(sce)$sample))
  outer_order <- c("PNF", "ANNUBP", "MPNST")

  expect_warning(
    g_empty <- plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "sample_type",
      group_order = sample_order,
      outer_order = outer_order,
      background = "sample_type",
      fg_subset_cells = quote(sample_type == "NOT_PRESENT")
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(g_empty, "patchwork")
  expect_equal(expected_nested_panel_count(sce, "sample", "sample_type"), 6L)
  expect_equal(expected_nested_ncol(sce, "sample", "sample_type"), 2L)

  g_partial <- plotOverlayPerNestedGroup(
    sce,
    foreground = "SOX10",
    group_col = "sample",
    outer_col = "sample_type",
    group_order = sample_order,
    outer_order = outer_order,
    background = "sample_type",
    fg_subset_cells = quote(sample_type == "MPNST")
  )
  expect_s3_class(g_partial, "patchwork")
  expect_equal(expected_nested_panel_count(sce, "sample", "sample_type"), 6L)
  expect_equal(expected_nested_ncol(sce, "sample", "sample_type"), 2L)
})

test_that("plotOverlayPerNestedGroup validates rasterization arguments early", {
  sce <- get_example_sce()

  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "sample_type",
      background = "solid",
      bg_raster = NA
    ),
    "`bg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "sample_type",
      background = "solid",
      fg_raster = c(TRUE, FALSE)
    ),
    "`fg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlayPerNestedGroup(
      sce,
      foreground = "SOX10",
      group_col = "sample",
      outer_col = "sample_type",
      background = "solid",
      raster_dpi = Inf
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
})
