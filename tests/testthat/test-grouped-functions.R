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

expected_nested_panel_count <- function(sce, group_col, outer_col) {
  group_vals <- as.character(SummarizedExperiment::colData(sce)[[group_col]])
  outer_vals <- as.character(SummarizedExperiment::colData(sce)[[outer_col]])
  valid <- !is.na(group_vals) & !is.na(outer_vals)
  groups_per_outer <- vapply(
    unique(outer_vals[valid]),
    function(outer) length(unique(group_vals[valid & outer_vals == outer])),
    integer(1)
  )
  return(length(groups_per_outer) * max(groups_per_outer))
}

expected_nested_ncol <- function(sce, group_col, outer_col) {
  group_vals <- as.character(SummarizedExperiment::colData(sce)[[group_col]])
  outer_vals <- as.character(SummarizedExperiment::colData(sce)[[outer_col]])
  valid <- !is.na(group_vals) & !is.na(outer_vals)
  groups_per_outer <- vapply(
    unique(outer_vals[valid]),
    function(outer) length(unique(group_vals[valid & outer_vals == outer])),
    integer(1)
  )
  return(max(groups_per_outer))
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
  g_bivar <- plotOverlayPerGroup(
    sce,
    foreground = c("numeric_fg", "numeric_y"),
    group_col = "sample_type",
    background = "solid",
    fg_type = "bivariate",
    fg_limits = "shared",
    fg_subset_cells = subset_idx
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
      fg_limits = c(0, 5)
    ),
    "Numeric foreground limits are not supported"
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
  expect_equal(
    patchwork_plot_count(g),
    expected_nested_panel_count(sce, "sample", "sample_type")
  )
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
  expect_equal(
    patchwork_plot_count(g_empty),
    expected_nested_panel_count(sce, "sample", "sample_type")
  )
  expect_equal(
    g_empty$patches$layout$ncol,
    expected_nested_ncol(sce, "sample", "sample_type")
  )

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
  expect_equal(
    patchwork_plot_count(g_partial),
    expected_nested_panel_count(sce, "sample", "sample_type")
  )
  expect_equal(
    g_partial$patches$layout$ncol,
    expected_nested_ncol(sce, "sample", "sample_type")
  )
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
