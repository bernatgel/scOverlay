get_foreground_layer_data <- function(p) {
  layer_data <- lapply(p$layers, function(layer) layer$data)
  foreground_layers <- Filter(
    function(x) is.data.frame(x) && "foreground" %in% colnames(x),
    layer_data
  )
  expect_true(length(foreground_layers) >= 1L)
  return(foreground_layers[[length(foreground_layers)]])
}

make_na_foreground_sce <- function() {
  assays <- list(
    logcounts = matrix(
      c(NA_real_, 0, 1, 2, NA_real_),
      nrow = 1L,
      dimnames = list("GeneA", paste0("cell", seq_len(5L)))
    )
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = data.frame(
      numeric_score = c(NA_real_, 0, 1, 2, NA_real_),
      categorical_state = c("A", "B", NA, "A", NA),
      cluster = c("1", "1", "2", "2", "3"),
      row.names = paste0("cell", seq_len(5L))
    )
  )
  SingleCellExperiment::reducedDim(sce, "TSNE") <- matrix(
    c(
      0, 0,
      1, 0,
      2, 0,
      3, 0,
      4, 0
    ),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(colnames(sce), c("TSNE1", "TSNE2"))
  )
  return(sce)
}

scale_names <- function(p) {
  return(vapply(p$scales$scales, `[[`, character(1), "name"))
}

test_that("plotOverlay returns ggplot for continuous and categorical foreground", {
  sce <- get_example_sce()

  p_minimal <- plotOverlay(
    sce = sce,
    foreground = "SOX10"
  )
  expect_s3_class(p_minimal, "ggplot")
  expect_equal(p_minimal$layers[[1]]$aes_params$colour, "gray70")

  p_cont <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10"
  )
  expect_s3_class(p_cont, "ggplot")

  p_cat <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "tumor_stage"
  )
  expect_s3_class(p_cat, "ggplot")
})

test_that("reduced dimension coordinates are validated explicitly", {
  sce <- get_example_sce()
  coords <- scOverlay:::.get_reduced_dim_coords(sce, "TSNE")

  expect_s3_class(coords, "data.frame")
  expect_named(coords, c("X", "Y"))
  expect_equal(nrow(coords), ncol(sce))

  expect_error(
    scOverlay:::.get_reduced_dim_coords(sce, "NOT_A_DIM"),
    "`reduced_dim` must be one of reducedDimNames\\(sce\\)"
  )

  SingleCellExperiment::reducedDim(sce, "ONE_COL") <- matrix(seq_len(ncol(sce)), ncol = 1L)
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "ONE_COL",
      background = "solid",
      foreground = "SOX10"
    ),
    "must contain at least two columns"
  )
})

test_that("fg_subset_cells filters foreground only and keeps full background", {
  sce <- get_example_sce()
  idx <- as.character(SummarizedExperiment::colData(sce)$sample_type) ==
    as.character(unique(SummarizedExperiment::colData(sce)$sample_type)[1])

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    fg_subset_cells = idx
  )

  expect_equal(nrow(p$data), ncol(sce))
  expect_equal(sum(!is.na(p$data$foreground)), sum(idx))
})

test_that(".eval_subset handles NULL and logical vectors", {
  sce <- get_example_sce()
  ncells <- ncol(sce)
  idx <- rep(c(TRUE, FALSE), length.out = ncells)

  expect_identical(scOverlay:::.eval_subset(NULL, sce), rep(TRUE, ncells))
  expect_identical(scOverlay:::.eval_subset(idx, sce), idx)

  idx_na <- idx
  idx_na[2] <- NA
  expected <- idx_na
  expected[is.na(expected)] <- FALSE
  expect_warning(
    result <- scOverlay:::.eval_subset(idx_na, sce),
    "contains NA values"
  )
  expect_identical(result, expected)

  expect_error(
    scOverlay:::.eval_subset(c(TRUE, FALSE), sce),
    "must have length ncol\\(sce\\)"
  )
})

test_that(".eval_subset handles character cell-name selections", {
  sce <- get_example_sce()
  colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  selected <- colnames(sce)[c(1, 3, 5)]

  expect_identical(
    scOverlay:::.eval_subset(selected, sce),
    colnames(sce) %in% selected
  )

  sce_no_names <- sce
  colnames(sce_no_names) <- NULL
  expect_error(
    scOverlay:::.eval_subset(selected, sce_no_names),
    "colnames\\(sce\\) is NULL"
  )

  expect_error(
    scOverlay:::.eval_subset(c(selected, "missing-cell"), sce),
    "missing-cell"
  )
})

test_that(".eval_subset validates function results", {
  sce <- get_example_sce()
  idx <- seq_len(ncol(sce)) <= 3

  expect_identical(
    scOverlay:::.eval_subset(function(sce) idx, sce),
    idx
  )
  expect_error(
    scOverlay:::.eval_subset(function(sce) c(TRUE, FALSE), sce),
    "must have length ncol\\(sce\\)"
  )
  expect_error(
    scOverlay:::.eval_subset(function(sce) seq_len(ncol(sce)), sce),
    "must return a logical vector"
  )
})

test_that(".eval_subset validates expression results", {
  sce <- get_example_sce()
  sample_type <- SummarizedExperiment::colData(sce)$sample_type
  expected <- sample_type == sample_type[1]

  expect_identical(
    scOverlay:::.eval_subset(quote(sample_type == sample_type[1]), sce),
    expected
  )
  expect_error(
    scOverlay:::.eval_subset(quote(c(TRUE, FALSE)), sce),
    "must have length ncol\\(sce\\)"
  )
  expect_error(
    scOverlay:::.eval_subset(quote(seq_len(ncol(sce))), sce),
    "could not be evaluated|must return a logical vector"
  )
  expect_error(
    scOverlay:::.eval_subset(quote(missing_col == 1), sce),
    "missing_col"
  )
})

test_that(".eval_subset warns when supported inputs match no cells", {
  sce <- get_example_sce()
  colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  none <- rep(FALSE, ncol(sce))

  expect_warning(
    scOverlay:::.eval_subset(none, sce),
    "did not match any cells"
  )
  expect_warning(
    scOverlay:::.eval_subset(character(0), sce),
    "did not match any cells"
  )
  expect_warning(
    scOverlay:::.eval_subset(function(sce) none, sce),
    "did not match any cells"
  )
  expect_warning(
    scOverlay:::.eval_subset(quote(rep(FALSE, length(sample_type))), sce),
    "did not match any cells"
  )
})

test_that(".eval_subset errors for unsupported input types", {
  sce <- get_example_sce()

  expect_error(
    scOverlay:::.eval_subset(1, sce),
    "must be NULL, a logical vector, a character vector of cell names, a function, or an expression"
  )
})

test_that(".order_foreground orders foreground rows", {
  coords_fg <- data.frame(
    id = seq_len(6),
    foreground = c(3, NA, 1, 2, NA, 3)
  )

  expect_identical(scOverlay:::.order_foreground(coords_fg, "input"), coords_fg)

  expect_identical(
    scOverlay:::.order_foreground(coords_fg, "ascending")$id,
    c(2L, 5L, 3L, 4L, 1L, 6L)
  )
  expect_identical(
    scOverlay:::.order_foreground(coords_fg, "descending")$id,
    c(2L, 5L, 1L, 6L, 4L, 3L)
  )

  set.seed(17)
  random_once <- scOverlay:::.order_foreground(coords_fg, "random")
  set.seed(17)
  random_twice <- scOverlay:::.order_foreground(coords_fg, "random")
  expect_identical(random_once, random_twice)
  expect_false(identical(random_once$id, coords_fg$id))
})

test_that(".order_foreground orders factors by factor levels", {
  coords_fg <- data.frame(
    id = seq_len(4),
    foreground = factor(
      c("medium", "low", "high", "medium"),
      levels = c("low", "medium", "high")
    )
  )

  expect_identical(
    scOverlay:::.order_foreground(coords_fg, "ascending")$id,
    c(2L, 1L, 4L, 3L)
  )
  expect_identical(
    scOverlay:::.order_foreground(coords_fg, "descending")$id,
    c(3L, 1L, 4L, 2L)
  )
})

test_that("bg_dimming controls white dimming layer opacity", {
  sce <- get_example_sce()

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    bg_dimming = 0
  )

  expect_equal(p$layers[[2]]$aes_params$alpha, 0)
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "cluster",
      foreground = "SOX10",
      bg_dimming = 1.5
    ),
    "bg_dimming"
  )
})

test_that("background and foreground legends can be controlled independently", {
  sce <- get_example_sce()

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    bg_legend = FALSE,
    fg_legend = TRUE,
    bg_legend_title = "Clusters",
    fg_legend_title = "Expression"
  )

  expect_equal(p$scales$scales[[1]]$name, "Clusters")
  expect_equal(p$scales$scales[[2]]$name, "Expression")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")

  p_no_fg <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "tumor_stage",
    bg_legend = TRUE,
    fg_legend = FALSE
  )

  expect_length(p_no_fg$scales$scales, 2)
  expect_s3_class(ggplot2::ggplot_build(p_no_fg), "ggplot_built")
})

test_that("solid background does not add a background legend scale", {
	sce <- get_example_sce()

	p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_legend = TRUE
  )

	expect_length(p$scales$scales, 1)
  expect_equal(p$scales$scales[[1]]$name, "SOX10")
})

test_that("solid background uses one fixed color from bg_palette", {
  sce <- get_example_sce()

  p_default <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10"
  )
  expect_equal(p_default$layers[[1]]$aes_params$colour, "gray70")

  p_custom <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_palette = "#111111"
  )
  expect_equal(p_custom$layers[[1]]$aes_params$colour, "#111111")

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      bg_palette = c("red", "blue")
    ),
    "`bg_palette`"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      bg_palette = "not_a_colour"
    ),
    "`bg_palette`"
  )
})

test_that("solid foreground draws selected cells with a fixed colour", {
  sce <- get_example_sce()
  colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  sample_type <- SummarizedExperiment::colData(sce)$sample_type
  selected <- sample_type == sample_type[1]
  selected_names <- colnames(sce)[selected]

  p_solid_solid <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "solid"
  )
  expect_s3_class(p_solid_solid, "ggplot")
  expect_equal(p_solid_solid$layers[[length(p_solid_solid$layers)]]$aes_params$colour, "red")
  expect_length(p_solid_solid$scales$scales, 0)
  expect_s3_class(ggplot2::ggplot_build(p_solid_solid), "ggplot_built")

  p_metadata_bg <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_subset_cells = quote(sample_type == sample_type[1]),
    fg_palette = "#E41A1C"
  )
  expect_s3_class(p_metadata_bg, "ggplot")
  expect_equal(p_metadata_bg$layers[[length(p_metadata_bg$layers)]]$aes_params$colour, "#E41A1C")
  expect_equal(nrow(p_metadata_bg$layers[[length(p_metadata_bg$layers)]]$data), sum(selected))
  expect_s3_class(ggplot2::ggplot_build(p_metadata_bg), "ggplot_built")

  p_gene_bg <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    bg_type = "continuous",
    foreground = "solid",
    fg_subset_cells = selected
  )
  expect_s3_class(p_gene_bg, "ggplot")
  expect_equal(nrow(p_gene_bg$layers[[length(p_gene_bg$layers)]]$data), sum(selected))
  expect_s3_class(ggplot2::ggplot_build(p_gene_bg), "ggplot_built")

  p_names <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_subset_cells = selected_names
  )
  expect_s3_class(p_names, "ggplot")
  expect_equal(nrow(p_names$layers[[length(p_names$layers)]]$data), length(selected_names))

  p_function <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_subset_cells = function(sce) SummarizedExperiment::colData(sce)$sample_type == sample_type[1]
  )
  expect_s3_class(p_function, "ggplot")
  expect_equal(nrow(p_function$layers[[length(p_function$layers)]]$data), sum(selected))
})

test_that("solid foreground validates colour and ignores legend settings", {
  sce <- get_example_sce()

  p_default <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_palette = NULL,
    fg_legend = TRUE
  )
  expect_equal(p_default$layers[[length(p_default$layers)]]$aes_params$colour, "red")
  expect_length(p_default$scales$scales, 1)

  p_no_legend <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_palette = "blue",
    fg_legend = FALSE
  )
  expect_equal(p_no_legend$layers[[length(p_no_legend$layers)]]$aes_params$colour, "blue")
  expect_length(p_no_legend$scales$scales, 1)

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "solid",
      fg_palette = c("red", "blue")
    ),
    "`fg_palette` must be NULL or a single valid colour when `foreground = \"solid\"`"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "solid",
      fg_palette = "not_a_colour"
    ),
    "`fg_palette` must be NULL or a single valid colour when `foreground = \"solid\"`"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "solid",
      fg_palette = 1
    ),
    "`fg_palette` must be NULL or a single valid colour when `foreground = \"solid\"`"
  )

  expect_equal(scOverlay:::.resolve_solid_fg_colour(NULL), "red")
  expect_equal(scOverlay:::.resolve_solid_fg_colour("#E41A1C"), "#E41A1C")
})

test_that("solid foreground supports rasterization and ordering controls", {
  sce <- get_example_sce()

  p_raster <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "solid",
    fg_raster = TRUE
  )
  expect_s3_class(p_raster, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p_raster), "ggplot_built")

  p_input <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "solid",
    fg_order = "input"
  )
  expect_s3_class(p_input, "ggplot")

  set.seed(17)
  p_random <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "solid",
    fg_order = "random"
  )
  expect_s3_class(p_random, "ggplot")
  expect_false(identical(
    p_random$layers[[length(p_random$layers)]]$data$X,
    p_input$layers[[length(p_input$layers)]]$data$X
  ))

  p_ascending <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "solid",
    fg_order = "ascending"
  )
  expect_s3_class(p_ascending, "ggplot")
  expect_identical(
    p_ascending$layers[[length(p_ascending$layers)]]$data$X,
    p_input$layers[[length(p_input$layers)]]$data$X
  )

  p_descending <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "solid",
    fg_order = "descending"
  )
  expect_s3_class(p_descending, "ggplot")
  expect_identical(
    p_descending$layers[[length(p_descending$layers)]]$data$X,
    p_input$layers[[length(p_input$layers)]]$data$X
  )
})

test_that("none layers return ggplots with stable empty point layers", {
  sce <- get_example_sce()
  sample_type <- SummarizedExperiment::colData(sce)$sample_type

  p_bg_none_gene <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "SOX10"
  )
  expect_s3_class(p_bg_none_gene, "ggplot")
  expect_length(p_bg_none_gene$layers, 3)
  expect_equal(nrow(p_bg_none_gene$layers[[1]]$data), 0)
  expect_length(p_bg_none_gene$scales$scales, 1)
  expect_equal(p_bg_none_gene$scales$scales[[1]]$name, "SOX10")
  expect_s3_class(ggplot2::ggplot_build(p_bg_none_gene), "ggplot_built")

  p_bg_none_categorical <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "cluster",
    fg_type = "categorical"
  )
  expect_s3_class(p_bg_none_categorical, "ggplot")
  expect_length(p_bg_none_categorical$layers, 3)
  expect_equal(nrow(p_bg_none_categorical$layers[[1]]$data), 0)
  expect_length(p_bg_none_categorical$scales$scales, 1)
  expect_s3_class(ggplot2::ggplot_build(p_bg_none_categorical), "ggplot_built")

  p_bg_none_solid <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "solid",
    fg_subset_cells = quote(sample_type == sample_type[1])
  )
  expect_s3_class(p_bg_none_solid, "ggplot")
  expect_length(p_bg_none_solid$layers, 3)
  expect_equal(nrow(p_bg_none_solid$layers[[1]]$data), 0)
  expect_equal(nrow(p_bg_none_solid$layers[[3]]$data), sum(sample_type == sample_type[1]))
  expect_s3_class(ggplot2::ggplot_build(p_bg_none_solid), "ggplot_built")

  p_fg_none_metadata <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "none"
  )
  expect_s3_class(p_fg_none_metadata, "ggplot")
  expect_length(p_fg_none_metadata$layers, 3)
  expect_equal(nrow(p_fg_none_metadata$layers[[3]]$data), 0)
  expect_length(p_fg_none_metadata$scales$scales, 1)
  expect_equal(p_fg_none_metadata$scales$scales[[1]]$name, "sample_type")

  p_fg_none_gene <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "none"
  )
  expect_s3_class(p_fg_none_gene, "ggplot")
  expect_length(p_fg_none_gene$layers, 3)
  expect_equal(nrow(p_fg_none_gene$layers[[3]]$data), 0)
  expect_length(p_fg_none_gene$scales$scales, 1)
  expect_equal(p_fg_none_gene$scales$scales[[1]]$name, "SOX10")

  p_fg_none_solid <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "none"
  )
  expect_s3_class(p_fg_none_solid, "ggplot")
  expect_length(p_fg_none_solid$layers, 3)
  expect_equal(nrow(p_fg_none_solid$layers[[3]]$data), 0)
  expect_length(p_fg_none_solid$scales$scales, 0)
  expect_s3_class(ggplot2::ggplot_build(p_fg_none_solid), "ggplot_built")

  p_both_none <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "none"
  )
  expect_s3_class(p_both_none, "ggplot")
  expect_length(p_both_none$layers, 3)
  expect_equal(nrow(p_both_none$layers[[1]]$data), 0)
  expect_equal(nrow(p_both_none$layers[[3]]$data), 0)
  expect_length(p_both_none$scales$scales, 0)
  expect_equal(p_both_none$labels$title, "")
  expect_s3_class(ggplot2::ggplot_build(p_both_none), "ggplot_built")

  p_both_none_title <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "none",
    title = "Empty panel"
  )
  expect_equal(p_both_none_title$labels$title, "Empty panel")
})

test_that("none layers ignore raster and scale settings safely", {
  sce <- get_example_sce()

  p_bg_raster <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "none",
    foreground = "SOX10",
    bg_raster = TRUE,
    bg_palette = "not_a_palette",
    bg_legend = FALSE
  )
  expect_s3_class(p_bg_raster, "ggplot")
  expect_length(p_bg_raster$scales$scales, 1)
  expect_s3_class(ggplot2::ggplot_build(p_bg_raster), "ggplot_built")

  p_fg_raster <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "none",
    fg_raster = TRUE,
    fg_palette = "not_a_palette",
    fg_legend = FALSE,
    fg_order = "descending",
    fg_subset_cells = quote(FALSE)
  )
  expect_s3_class(p_fg_raster, "ggplot")
  expect_length(p_fg_raster$scales$scales, 1)
  expect_s3_class(ggplot2::ggplot_build(p_fg_raster), "ggplot_built")
})

test_that("empty foreground subsets still return plots", {
  sce <- get_example_sce()
  empty_subset <- rep(FALSE, ncol(sce))

  expect_warning(
    p_gene <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "sample_type",
      foreground = "SOX10",
      fg_subset_cells = empty_subset
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(p_gene, "ggplot")
  expect_length(p_gene$layers, 3)
  expect_equal(nrow(p_gene$layers[[3]]$data), 0)
  expect_s3_class(ggplot2::ggplot_build(p_gene), "ggplot_built")

  expect_warning(
    p_categorical <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "sample_type",
      foreground = "cluster",
      fg_type = "categorical",
      fg_subset_cells = quote(sample_type == "NOT_PRESENT")
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(p_categorical, "ggplot")
  expect_equal(nrow(p_categorical$layers[[3]]$data), 0)
  expect_s3_class(ggplot2::ggplot_build(p_categorical), "ggplot_built")

  expect_warning(
    p_solid <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "solid",
      fg_subset_cells = empty_subset
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(p_solid, "ggplot")
  expect_equal(nrow(p_solid$layers[[3]]$data), 0)

  expect_warning(
    p_bg_none <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "none",
      foreground = "SOX10",
      fg_subset_cells = empty_subset
    ),
    "drawing an empty foreground layer"
  )
  expect_s3_class(p_bg_none, "ggplot")
  expect_equal(nrow(p_bg_none$layers[[3]]$data), 0)
  expect_s3_class(ggplot2::ggplot_build(p_bg_none), "ggplot_built")

  p_fg_none <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "sample_type",
    foreground = "none",
    fg_subset_cells = empty_subset
  )
  expect_s3_class(p_fg_none, "ggplot")
  expect_equal(nrow(p_fg_none$layers[[3]]$data), 0)
  expect_s3_class(ggplot2::ggplot_build(p_fg_none), "ggplot_built")
})

test_that("background variables support categorical and continuous scales", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_bg <- rep(c(1, 2), length.out = ncol(sce))

  p_categorical_coldata <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10"
  )
  expect_s3_class(p_categorical_coldata, "ggplot")
  expect_equal(p_categorical_coldata$scales$scales[[1]]$name, "cluster")
  expect_s3_class(ggplot2::ggplot_build(p_categorical_coldata), "ggplot_built")

  p_numeric_coldata <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "numeric_bg",
    foreground = "SOX10"
  )
  expect_s3_class(p_numeric_coldata, "ggplot")
  expect_equal(p_numeric_coldata$scales$scales[[1]]$name, "numeric_bg")
  expect_s3_class(ggplot2::ggplot_build(p_numeric_coldata), "ggplot_built")

  p_feature_auto <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "S100B"
  )
  expect_s3_class(p_feature_auto, "ggplot")
  expect_equal(p_feature_auto$scales$scales[[1]]$name, "SOX10")
  expect_s3_class(ggplot2::ggplot_build(p_feature_auto), "ggplot_built")

  p_feature_continuous <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "S100B",
    bg_type = "continuous"
  )
  expect_s3_class(p_feature_continuous, "ggplot")
  expect_equal(p_feature_continuous$scales$scales[[1]]$name, "SOX10")
  expect_s3_class(ggplot2::ggplot_build(p_feature_continuous), "ggplot_built")

  p_feature_categorical <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "S100B",
    bg_type = "categorical",
    bg_palette = function(n) grDevices::hcl.colors(n, "Dark 3")
  )
  expect_s3_class(p_feature_categorical, "ggplot")
  expect_equal(p_feature_categorical$scales$scales[[1]]$name, "SOX10")
  expect_s3_class(ggplot2::ggplot_build(p_feature_categorical), "ggplot_built")

  p_numeric_categorical <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "numeric_bg",
    foreground = "SOX10",
    bg_type = "categorical"
  )
  expect_s3_class(p_numeric_categorical, "ggplot")
  expect_equal(p_numeric_categorical$scales$scales[[1]]$name, "numeric_bg")
  expect_s3_class(ggplot2::ggplot_build(p_numeric_categorical), "ggplot_built")
})

test_that("layer colour scale limits are validated", {
  expect_silent(scOverlay:::.validate_layer_limits(NULL, "fg_limits"))
  expect_silent(scOverlay:::.validate_layer_limits(c(0, 1), "fg_limits"))
  expect_silent(scOverlay:::.validate_layer_limits("shared", "fg_limits"))
  expect_silent(scOverlay:::.validate_layer_limits(NULL, "bg_limits"))
  expect_silent(scOverlay:::.validate_layer_limits(c(0, 1), "bg_limits"))
  expect_silent(scOverlay:::.validate_layer_limits("shared", "bg_limits"))

  expect_error(
    scOverlay:::.validate_layer_limits("global", "fg_limits"),
    "`fg_limits` must be NULL, \"shared\", or a numeric vector of length 2"
  )
  expect_error(
    scOverlay:::.validate_layer_limits(c(0, 1, 2), "bg_limits"),
    "`bg_limits` must be NULL, \"shared\", or a numeric vector of length 2"
  )
  expect_error(
    scOverlay:::.validate_layer_limits(c(0, NA), "fg_limits"),
    "`fg_limits` must contain two finite non-NA values"
  )
  expect_error(
    scOverlay:::.validate_layer_limits(c(0, Inf), "bg_limits"),
    "`bg_limits` must contain two finite non-NA values"
  )
  expect_error(
    scOverlay:::.validate_layer_limits(c(1, 1), "fg_limits"),
    "`fg_limits` must be strictly increasing"
  )
  expect_error(
    scOverlay:::.validate_layer_limits(c(2, 1), "bg_limits"),
    "`bg_limits` must be strictly increasing"
  )
})

test_that("plotOverlay supports foreground and background limits", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_bg <- seq_len(ncol(sce))

  p_fg_shared <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    fg_limits = "shared"
  )
  expect_s3_class(p_fg_shared, "ggplot")
  expect_length(p_fg_shared$scales$scales[[1]]$limits, 2)
  expect_s3_class(ggplot2::ggplot_build(p_fg_shared), "ggplot_built")

  p_bg_shared <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "numeric_bg",
    foreground = "SOX10",
    bg_limits = "shared"
  )
  expect_s3_class(p_bg_shared, "ggplot")
  expect_equal(p_bg_shared$scales$scales[[1]]$limits, c(1, ncol(sce)))
  expect_s3_class(ggplot2::ggplot_build(p_bg_shared), "ggplot_built")

  p_bg_numeric <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "numeric_bg",
    foreground = "SOX10",
    bg_limits = c(0, 20)
  )
  expect_equal(p_bg_numeric$scales$scales[[1]]$limits, c(0, 20))
  expect_s3_class(ggplot2::ggplot_build(p_bg_numeric), "ggplot_built")
})

test_that("limits supplied to non-continuous layers warn and are ignored", {
  sce <- get_example_sce()

  expect_warning(
    p_fg_categorical <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "cluster",
      fg_type = "categorical",
      fg_limits = c(0, 1)
    ),
    "`fg_limits` is ignored"
  )
  expect_s3_class(p_fg_categorical, "ggplot")

  expect_warning(
    p_fg_solid <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "solid",
      fg_limits = "shared"
    ),
    "`fg_limits` is ignored"
  )
  expect_s3_class(p_fg_solid, "ggplot")

  expect_warning(
    p_fg_none <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "none",
      fg_limits = c(0, 1)
    ),
    "`fg_limits` is ignored"
  )
  expect_s3_class(p_fg_none, "ggplot")

  expect_warning(
    p_bg_categorical <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "cluster",
      foreground = "SOX10",
      bg_limits = c(0, 1)
    ),
    "`bg_limits` is ignored"
  )
  expect_s3_class(p_bg_categorical, "ggplot")

  expect_warning(
    p_bg_solid <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      bg_limits = "shared"
    ),
    "`bg_limits` is ignored"
  )
  expect_s3_class(p_bg_solid, "ggplot")

  expect_warning(
    p_bg_none <- plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "none",
      foreground = "SOX10",
      bg_limits = c(0, 1)
    ),
    "`bg_limits` is ignored"
  )
  expect_s3_class(p_bg_none, "ggplot")
})

test_that("background validation is explicit", {
  sce <- get_example_sce()

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "cluster",
      foreground = "SOX10",
      bg_type = "continuous"
    ),
    "bg_type = 'continuous'.*numeric background"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "missing_background",
      foreground = "SOX10"
    ),
    "`background` must be \"solid\", \"none\", a feature in rownames\\(sce\\), or a column in colData\\(sce\\)"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "SOX10",
      foreground = "S100B",
      bg_assay = "missing_assay"
    ),
    "`bg_assay` must be one of assayNames\\(sce\\)"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "SOX10",
      foreground = "S100B",
      bg_assay = 1
    ),
    "`bg_assay` must be a single non-NA character string"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "SOX10",
      foreground = "S100B",
      bg_assay = c("counts", "logcounts")
    ),
    "`bg_assay` must be a single non-NA character string"
  )

  p_coldata_ignores_bg_assay <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    bg_assay = "missing_assay"
  )
  expect_s3_class(p_coldata_ignores_bg_assay, "ggplot")
})

test_that("background and foreground scales remain independent", {
  sce <- get_example_sce()

  p_gene_gene <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "S100B"
  )
  expect_s3_class(p_gene_gene, "ggplot")
  expect_equal(vapply(p_gene_gene$scales$scales, `[[`, character(1), "name"), c("SOX10", "S100B"))
  expect_s3_class(ggplot2::ggplot_build(p_gene_gene), "ggplot_built")

  p_gene_categorical <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "cluster",
    fg_type = "categorical"
  )
  expect_s3_class(p_gene_categorical, "ggplot")
  expect_equal(vapply(p_gene_categorical$scales$scales, `[[`, character(1), "name"), c("SOX10", "cluster"))
  expect_s3_class(ggplot2::ggplot_build(p_gene_categorical), "ggplot_built")
})

test_that("feature backgrounds support rasterization", {
  sce <- get_example_sce()

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "SOX10",
    foreground = "S100B",
    bg_raster = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that(".geom_point_overlay returns vector and raster point layers", {
  layer_vector <- scOverlay:::.geom_point_overlay(
    data = data.frame(x = 1, y = 1),
    mapping = ggplot2::aes(x, y),
    raster = FALSE
  )
  layer_raster <- scOverlay:::.geom_point_overlay(
    data = data.frame(x = 1, y = 1),
    mapping = ggplot2::aes(x, y),
    raster = TRUE,
    raster_dpi = 600
  )

  expect_s3_class(layer_vector, "LayerInstance")
  expect_s3_class(layer_raster, "LayerInstance")
  expect_true(is.function(layer_vector$draw_geom))
  expect_true(is.function(layer_raster$draw_geom))
})

test_that("plotOverlay supports optional rasterized point layers", {
  sce <- get_example_sce()

  p_vector <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_raster = FALSE,
    fg_raster = FALSE
  )
  expect_s3_class(p_vector, "ggplot")

  p_bg <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_raster = TRUE,
    fg_raster = FALSE
  )
  expect_s3_class(p_bg, "ggplot")

  p_fg <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_raster = FALSE,
    fg_raster = TRUE
  )
  expect_s3_class(p_fg, "ggplot")

  p_both <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_raster = TRUE,
    fg_raster = TRUE
  )
  expect_s3_class(p_both, "ggplot")

  p_dpi <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    bg_raster = TRUE,
    fg_raster = TRUE,
    raster_dpi = 600
  )
  expect_s3_class(p_dpi, "ggplot")
})

test_that("plotOverlay validates rasterization arguments", {
  sce <- get_example_sce()

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      bg_raster = NA
    ),
    "`bg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      fg_raster = c(TRUE, FALSE)
    ),
    "`fg_raster` must be TRUE or FALSE"
  )
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      raster_dpi = 0
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      raster_dpi = NA_real_
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      raster_dpi = Inf
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "SOX10",
      raster_dpi = -Inf
    ),
    "`raster_dpi` must be a single positive finite numeric value"
  )
})

test_that("fg_order works through plotOverlay", {
	sce <- get_example_sce()

  p_ascending <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    fg_order = "ascending",
    fg_legend = TRUE
  )
  expect_s3_class(p_ascending, "ggplot")
  expect_length(p_ascending$scales$scales, 1)
  expect_equal(p_ascending$scales$scales[[1]]$name, "SOX10")
  expect_s3_class(ggplot2::ggplot_build(p_ascending), "ggplot_built")

  set.seed(17)
  p_random <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    fg_order = "random"
  )
  expect_s3_class(p_random, "ggplot")
})

test_that("NA foreground values are dropped from the foreground layer", {
  sce <- get_example_sce()
  fg <- seq_len(ncol(sce))
  fg[c(2, 5, 8)] <- NA_integer_
  SummarizedExperiment::colData(sce)$fg_with_na <- fg

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "fg_with_na",
    fg_order = "input"
  )

  foreground_layer <- p$layers[[length(p$layers)]]$data
  expect_equal(nrow(foreground_layer), sum(!is.na(fg)))
  expect_false(any(is.na(foreground_layer$foreground)))
})

test_that("NA foreground values are dropped before ordering", {
  sce <- make_na_foreground_sce()

  p_numeric_ascending <- plotOverlay(
    sce = sce,
    foreground = "numeric_score",
    background = "solid",
    reduced_dim = "TSNE",
    fg_type = "continuous",
    fg_order = "ascending"
  )
  numeric_ascending <- get_foreground_layer_data(p_numeric_ascending)
  expect_equal(nrow(numeric_ascending), 3L)
  expect_false(any(is.na(numeric_ascending$foreground)))
  expect_equal(numeric_ascending$foreground, c(0, 1, 2))

  p_numeric_descending <- plotOverlay(
    sce = sce,
    foreground = "numeric_score",
    background = "solid",
    reduced_dim = "TSNE",
    fg_type = "continuous",
    fg_order = "descending"
  )
  numeric_descending <- get_foreground_layer_data(p_numeric_descending)
  expect_equal(nrow(numeric_descending), 3L)
  expect_false(any(is.na(numeric_descending$foreground)))
  expect_equal(numeric_descending$foreground, c(2, 1, 0))

  p_gene <- plotOverlay(
    sce = sce,
    foreground = "GeneA",
    background = "solid",
    reduced_dim = "TSNE",
    fg_order = "ascending"
  )
  gene_foreground <- get_foreground_layer_data(p_gene)
  expect_equal(nrow(gene_foreground), 3L)
  expect_false(any(is.na(gene_foreground$foreground)))
  expect_equal(gene_foreground$foreground, c(0, 1, 2))

  p_categorical <- plotOverlay(
    sce = sce,
    foreground = "categorical_state",
    background = "solid",
    reduced_dim = "TSNE",
    fg_type = "categorical"
  )
  categorical_foreground <- get_foreground_layer_data(p_categorical)
  expect_equal(nrow(categorical_foreground), 3L)
  expect_false(any(is.na(categorical_foreground$foreground)))
  expect_equal(categorical_foreground$foreground, c("A", "B", "A"))
})

test_that("solid and none foregrounds bypass variable NA filtering", {
  sce <- make_na_foreground_sce()

  p_solid <- plotOverlay(
    sce = sce,
    foreground = "solid",
    background = "solid",
    reduced_dim = "TSNE"
  )
  solid_foreground <- get_foreground_layer_data(p_solid)
  expect_s3_class(p_solid, "ggplot")
  expect_equal(nrow(solid_foreground), ncol(sce))
  expect_equal(solid_foreground$foreground, rep("all", ncol(sce)))

  p_none <- plotOverlay(
    sce = sce,
    foreground = "none",
    background = "solid",
    reduced_dim = "TSNE"
  )
  none_foreground <- get_foreground_layer_data(p_none)
  expect_s3_class(p_none, "ggplot")
  expect_equal(nrow(none_foreground), 0L)
})

test_that("empty foreground subsets still produce empty foreground layers with NA-capable data", {
  sce <- make_na_foreground_sce()

  expect_warning(
    p <- plotOverlay(
      sce = sce,
      foreground = "numeric_score",
      background = "solid",
      reduced_dim = "TSNE",
      fg_subset_cells = rep(FALSE, ncol(sce))
    ),
    "drawing an empty foreground layer"
  )
  foreground <- get_foreground_layer_data(p)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(foreground), 0L)
})

test_that("fg_type continuous requires numeric foreground data", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_fg <- seq_len(ncol(sce))
  SummarizedExperiment::colData(sce)$character_fg <- rep(c("A", "B"), length.out = ncol(sce))
  SummarizedExperiment::colData(sce)$factor_fg <- factor(
    rep(c("low", "high"), length.out = ncol(sce)),
    levels = c("low", "high")
  )

  p_numeric <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "numeric_fg",
    fg_type = "continuous"
  )
  expect_s3_class(p_numeric, "ggplot")

  p_gene <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "SOX10",
    fg_type = "continuous"
  )
  expect_s3_class(p_gene, "ggplot")

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "character_fg",
      fg_type = "continuous"
    ),
    "fg_type = 'continuous'.*numeric foreground"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "factor_fg",
      fg_type = "continuous"
    ),
    "fg_type = 'continuous'.*numeric foreground"
  )
})

test_that("fg_type categorical remains permissive", {
  sce <- get_example_sce()
  SummarizedExperiment::colData(sce)$numeric_fg <- rep(c(1, 2), length.out = ncol(sce))
  SummarizedExperiment::colData(sce)$character_fg <- rep(c("A", "B"), length.out = ncol(sce))
  SummarizedExperiment::colData(sce)$factor_fg <- factor(
    rep(c("low", "high"), length.out = ncol(sce)),
    levels = c("low", "high")
  )

  p_character <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "character_fg",
    fg_type = "categorical"
  )
  expect_s3_class(p_character, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p_character), "ggplot_built")

  p_factor <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "factor_fg",
    fg_type = "categorical"
  )
  expect_s3_class(p_factor, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p_factor), "ggplot_built")

  p_numeric <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = "numeric_fg",
    fg_type = "categorical"
  )
  expect_s3_class(p_numeric, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p_numeric), "ggplot_built")

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "character_fg",
      fg_type = "discrete"
    ),
    "'arg' should be one of"
  )
})

test_that("fg_type bivariate uses two numeric foreground variables", {
  sce <- make_na_foreground_sce()
  SummarizedExperiment::colData(sce)$score_y <- c(0, 1, NA, 3, 4)

  p <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = c("numeric_score", "score_y"),
    fg_type = "bivariate",
    fg_limits = NULL
  )
  foreground <- get_foreground_layer_data(p)
  expect_s3_class(p, "ggplot")
  expect_true(".fg_colour" %in% colnames(foreground))
  expect_equal(nrow(foreground), 2L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", foreground$.fg_colour)))

  p_shared <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "solid",
    foreground = c("numeric_score", "score_y"),
    fg_type = "bivariate",
    fg_limits = "shared"
  )
  expect_s3_class(p_shared, "ggplot")

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = c("numeric_score", "score_y"),
      fg_type = "bivariate",
      fg_limits = c(0, 5)
    ),
    "Numeric foreground limits are not supported"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = c("numeric_score", "score_y")
    ),
    "`foreground` must be a single non-NA character string"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = "numeric_score",
      fg_type = "bivariate"
    ),
    "character vector of length 2"
  )

  expect_error(
    plotOverlay(
      sce,
      reduced_dim = "TSNE",
      background = "solid",
      foreground = c("numeric_score", "categorical_state"),
      fg_type = "bivariate"
    ),
    "requires both foreground variables to be numeric"
  )
})

test_that("plotOverlay validates key inputs", {
  sce <- get_example_sce()

  expect_error(
    plotOverlay(sce = data.frame(), reduced_dim = "TSNE", background = "cluster", foreground = "SOX10"),
    "`sce` must be a SingleCellExperiment object"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "NOT_A_DIM", background = "cluster", foreground = "SOX10"),
    "`reduced_dim` must be one of reducedDimNames\\(sce\\)"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "NOT_A_COL", foreground = "SOX10"),
    "`background`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = 1, foreground = "SOX10"),
    "`background`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = c("cluster", "sample"), foreground = "SOX10"),
    "`background`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "NOT_A_FEATURE"),
    "`foreground`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = 1),
    "`foreground`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = c("SOX10", "S100B")),
    "`foreground`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", bg_legend = NA),
    "`bg_legend`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", fg_legend = NA),
    "`fg_legend`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "gray", foreground = "SOX10"),
    "`background`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", bg_dimming = Inf),
    "`bg_dimming`"
  )

  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", fg_assay = 1),
    "`fg_assay`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", fg_assay = c("counts", "logcounts")),
    "`fg_assay`"
  )
  expect_error(
    plotOverlay(sce = sce, reduced_dim = "TSNE", background = "cluster", foreground = "SOX10", fg_assay = "missing_assay"),
    "`fg_assay`"
  )
})
