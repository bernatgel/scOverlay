expect_valid_colours <- function(cols) {
  valid <- vapply(cols, function(col) {
    !inherits(try(grDevices::col2rgb(col), silent = TRUE), "try-error")
  }, logical(1))
  expect_true(all(valid))
  return(invisible(TRUE))
}

test_that("getPalette resolves scOverlay defaults", {
  expect_false("default" %in% names(scOverlay:::scOverlay_palettes$categorical))
  expect_false("default" %in% names(scOverlay:::scOverlay_palettes$continuous))

  pal_continuous <- getPalette(NULL, type = "continuous", n = 10)
  expect_type(pal_continuous, "character")
  expect_length(pal_continuous, 10)
  expect_equal(pal_continuous, getPalette("gray_red", type = "continuous", n = 10))
  expect_valid_colours(pal_continuous)

  pal_categorical <- getPalette(NULL, type = "categorical", n = 5)
  expect_type(pal_categorical, "character")
  expect_length(pal_categorical, 5)
  expect_equal(pal_categorical, unname(scOverlay:::scOverlay_palettes$categorical$scOverlay[seq_len(5)]))
  expect_valid_colours(pal_categorical)
})

test_that("mode-specific palette defaults remain unchanged", {
  expect_identical(scOverlay:::.default_palette_name("categorical"), "scOverlay")
  expect_identical(scOverlay:::.default_palette_name("continuous"), "gray_red")
  expect_identical(scOverlay:::.default_palette_name("bivariate"), "bivariateBright")
  expect_identical(scOverlay:::.default_palette_name("combination"), "combinationBright")

  expect_identical(
    scOverlay:::.resolve_fg_palette(NULL, "bivariate"),
    "bivariateBright"
  )
  expect_identical(
    scOverlay:::.resolve_fg_palette(NULL, "combination"),
    "combinationBright"
  )
  expect_null(scOverlay:::.resolve_fg_palette(NULL, "continuous"))
})

test_that("bivariate palette corners and orientation remain unchanged", {
  expected <- list(
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
  )

  expect_identical(scOverlay:::scOverlay_palettes$bivariate, expected)
  expect_identical(
    scOverlay:::.bivariate_palette_corners(),
    expected$bivariateBright
  )

  for (palette_name in names(expected)) {
    corners <- expected[[palette_name]]
    expect_identical(
      scOverlay:::.bivariate_palette_corners(palette_name),
      corners
    )

    palette_matrix <- scOverlay:::.bivariate_palette(palette_name, n = 5L)
    matrix_corners <- c(
      low_xy = palette_matrix[1L, 1L],
      high_x = palette_matrix[1L, 5L],
      high_y = palette_matrix[5L, 1L],
      high_xy = palette_matrix[5L, 5L]
    )
    expect_equal(
      unname(grDevices::col2rgb(matrix_corners)),
      unname(grDevices::col2rgb(corners))
    )
  }
})

test_that("combinationBright ordering remains unchanged", {
  expected <- c(
    "#D9D9D9", "#ff4346", "#55b2ff", "#68d865",
    "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
    "#F781BF", "#00CED1", "#B3DE69", "#FDB462",
    "#80B1D3", "#FB8072", "#BC80BD", "#8DD3C7",
    "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02",
    "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666", "#1F78B4", "#33A02C",
    "#FB9A99", "#CAB2D6", "#FDBF6F", "#B2DF8A"
  )

  expect_identical(
    scOverlay:::scOverlay_palettes$categorical$combinationBright,
    expected
  )
  expect_identical(
    getPalette("combinationBright", type = "categorical", n = 32L),
    expected
  )
})

test_that("getPalette resolves named scOverlay palettes by type", {
  pal_gray_red <- getPalette("gray_red", type = "continuous", n = 256)
  expect_type(pal_gray_red, "character")
  expect_length(pal_gray_red, 256)
  expect_valid_colours(pal_gray_red)

  pal_gray_red_short <- getPalette("gray_red", type = "continuous", n = 5)
  expect_type(pal_gray_red_short, "character")
  expect_length(pal_gray_red_short, 5)
  expect_valid_colours(pal_gray_red_short)

  pal_soft <- getPalette("soft", type = "categorical", n = 10)
  expect_type(pal_soft, "character")
  expect_length(pal_soft, 10)
  expect_valid_colours(pal_soft)

  pal_combination_15 <- getPalette("combinationBright", type = "categorical", n = 15)
  expect_type(pal_combination_15, "character")
  expect_length(pal_combination_15, 15)
  expect_equal(pal_combination_15[[1]], "#D9D9D9")
  expect_true(all(grepl("^#", pal_combination_15)))
  expect_valid_colours(pal_combination_15)

  pal_combination_32 <- getPalette("combinationBright", type = "categorical", n = 32)
  expect_length(pal_combination_32, 32)
  expect_equal(pal_combination_32[[1]], "#D9D9D9")
  expect_equal(pal_combination_32[[2]], "#ff4346")
  expect_equal(pal_combination_32[[3]], "#55b2ff")
  expect_equal(pal_combination_32[[4]], "#68d865")
  expect_valid_colours(pal_combination_32)
})

test_that("getPalette interpolates user-supplied continuous color vectors", {
  pal_vec <- getPalette(c("#E6E6E6", "#67000D"), type = "continuous", n = 10)
  expect_type(pal_vec, "character")
  expect_length(pal_vec, 10)
  expect_valid_colours(pal_vec)

  pal_non_integer <- getPalette("gray_red", type = "continuous", n = 5.8)
  expect_length(pal_non_integer, 5)

  pal_categorical_non_integer <- getPalette("scOverlay", type = "categorical", n = 5.8)
  expect_length(pal_categorical_non_integer, 5)
})

test_that("getPalette resolves external palette names", {
  pal_brewer <- getPalette("Set1", type = "categorical", n = 3, values = c("A", "B", "C"))
  expect_type(pal_brewer, "character")
  expect_length(pal_brewer, 3)
  expect_equal(names(pal_brewer), c("A", "B", "C"))
  expect_valid_colours(pal_brewer)

  viridis_names <- c(
    "viridis",
    "magma",
    "inferno",
    "plasma",
    "cividis",
    "rocket",
    "mako",
    "turbo"
  )
  for (pal in viridis_names) {
    cols <- getPalette(pal, type = "continuous", n = 5)
    expect_type(cols, "character")
    expect_length(cols, 5)
    expect_true(all(grepl("^#", cols)))
    expect_valid_colours(cols)
  }
})

test_that("ColorBrewer palettes handle supported and out-of-range sizes", {
  brewer_palettes <- c(
    Blues = "seq",
    Set1 = "qual",
    RdBu = "div"
  )

  for (palette_name in names(brewer_palettes)) {
    palette_info <- RColorBrewer::brewer.pal.info[palette_name, , drop = FALSE]
    maximum_colours <- as.integer(palette_info$maxcolors[[1L]])
    expect_identical(
      as.character(palette_info$category[[1L]]),
      unname(brewer_palettes[[palette_name]])
    )

    for (n in unique(c(1L, 2L, 3L, maximum_colours))) {
      expected <- RColorBrewer::brewer.pal(max(3L, n), palette_name)
      if (n < 3L && identical(brewer_palettes[[palette_name]], "qual")) {
        expected <- expected[seq_len(n)]
      } else if (n == 1L) {
        expected <- expected[[2L]]
      } else if (n == 2L) {
        expected <- expected[c(1L, 3L)]
      }
      expect_silent(
        cols <- scOverlay:::.brewer_palette(
          palette_name,
          type = "continuous",
          n = n
        )
      )
      expect_equal(cols, expected)
      expect_length(cols, n)
      expect_valid_colours(cols)
    }

    expect_silent(
      base_cols <- scOverlay:::.brewer_palette(
        palette_name,
        type = "categorical",
        n = maximum_colours + 3L
      )
    )
    expect_equal(
      base_cols,
      RColorBrewer::brewer.pal(maximum_colours, palette_name)
    )
    expect_length(base_cols, maximum_colours)
    expect_valid_colours(base_cols)
  }
})

test_that("getPalette uses deliberate small ColorBrewer selections", {
  brewer_palettes <- c(
    Blues = "seq",
    Set1 = "qual",
    RdBu = "div"
  )

  for (palette_name in names(brewer_palettes)) {
    base_cols <- RColorBrewer::brewer.pal(3L, palette_name)
    if (identical(brewer_palettes[[palette_name]], "qual")) {
      expected_one <- base_cols[[1L]]
      expected_two <- base_cols[c(1L, 2L)]
    } else {
      expected_one <- base_cols[[2L]]
      expected_two <- base_cols[c(1L, 3L)]
    }

    for (type in c("continuous", "categorical")) {
      expect_silent(
        one_col <- getPalette(palette_name, type = type, n = 1L)
      )
      expect_equal(one_col, expected_one)
      expect_valid_colours(one_col)

      expect_silent(
        two_cols <- getPalette(palette_name, type = type, n = 2L)
      )
      expect_equal(two_cols, expected_two)
      expect_valid_colours(two_cols)
    }
  }

  expect_equal(
    getPalette("RdBu", type = "continuous", n = 2L),
    RColorBrewer::brewer.pal(3L, "RdBu")[c(1L, 3L)]
  )
})

test_that("getPalette finalizes ColorBrewer palettes for either requested type", {
  brewer_palettes <- c("Blues", "Set1", "RdBu")

  for (palette_name in brewer_palettes) {
    maximum_colours <- as.integer(
      RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    )
    oversized_n <- maximum_colours + 3L

    expect_silent(
      continuous_cols <- getPalette(
        palette_name,
        type = "continuous",
        n = oversized_n
      )
    )
    expect_length(continuous_cols, oversized_n)
    expect_valid_colours(continuous_cols)

    expect_silent(
      categorical_cols <- getPalette(
        palette_name,
        type = "categorical",
        n = maximum_colours
      )
    )
    expect_length(categorical_cols, maximum_colours)
    expect_valid_colours(categorical_cols)

    expect_warning(
      oversized_categorical_cols <- getPalette(
        palette_name,
        type = "categorical",
        n = oversized_n
      ),
      "Categorical palette has fewer colours than requested"
    )
    expect_length(oversized_categorical_cols, oversized_n)
    expect_valid_colours(oversized_categorical_cols)
  }
})

test_that("getPalette honors named categorical colors when values match names", {
  values <- c("3", "1", "2")
  custom <- c("1" = "#05d165", "2" = "#e6194b", "3" = "#ebad00")

  pal <- getPalette(
    palette = custom,
    type = "categorical",
    n = 3,
    values = values
  )

  expect_equal(names(pal), values)
  expect_equal(unname(pal), c("#ebad00", "#05d165", "#e6194b"))
})

test_that("getPalette handles named categorical palette matching strictly", {
  expect_equal(
    getPalette(c(A = "red", B = "blue"), type = "categorical", values = c("B", "A")),
    c(B = "blue", A = "red")
  )

  expect_equal(
    getPalette(c("1" = "red", "2" = "blue"), type = "categorical", values = c("A", "B")),
    c(A = "red", B = "blue")
  )

  expect_error(
    getPalette(c(A = "red", B = "blue"), type = "categorical", values = c("A", "C")),
    "does not define colours for: C",
    fixed = TRUE
  )

  expect_error(
    getPalette("scOverlay", type = "categorical", values = c("1", "2", "Schwann")),
    "Schwann"
  )

  cluster_values <- c("2", "1")
  cluster_pal <- getPalette("scOverlay", type = "categorical", values = cluster_values)
  expect_identical(names(cluster_pal), cluster_values)
  expect_equal(
    unname(cluster_pal),
    unname(scOverlay:::scOverlay_palettes$categorical$scOverlay[cluster_values])
  )

  biological_values <- c("Schwann", "Fibroblast")
  biological_pal <- getPalette("scOverlay", type = "categorical", values = biological_values)
  expect_identical(names(biological_pal), biological_values)
  expect_equal(length(biological_pal), length(biological_values))
  expect_equal(
    unname(biological_pal),
    unname(scOverlay:::scOverlay_palettes$categorical$scOverlay[seq_along(biological_values)])
  )
})

test_that("getPalette warns and interpolates short categorical palettes", {
  expect_warning(
    pal <- getPalette(c("#000000", "#ffffff"), type = "categorical", n = 4),
    "fewer colours"
  )
  expect_type(pal, "character")
  expect_length(pal, 4)
  expect_valid_colours(pal)
})

test_that("getPalette errors for unknown names and invalid colors", {
  expect_error(getPalette("not_a_palette", type = "continuous", n = 5), "Unknown palette")
  expect_error(getPalette(c("red", "not_a_colour"), type = "categorical", n = 2), "invalid colours")
  expect_error(getPalette(1, type = "continuous", n = 5), "`palette`")
  expect_error(getPalette(character(0), type = "continuous", n = 5), "resolve")
  expect_error(getPalette(c("not_a_colour"), type = "categorical", n = 1), "Unknown palette")
  expect_error(getPalette("not_a_palette", type = "continuous", n = 0), "`n`")
  expect_error(getPalette(type = "discrete"), "'arg' should be one of")
})

test_that("continuous palettes are unaffected by categorical name matching", {
  pal <- getPalette(c(A = "gray90", B = "red"), type = "continuous", n = 5, values = c("A", "C"))
  expect_type(pal, "character")
  expect_length(pal, 5)
  expect_null(names(pal))
  expect_valid_colours(pal)
})

test_that("listPalettes returns palette metadata", {
  expected_bivariate <- c(
    "bivariateBright",
    "bivariateOrangeBlue",
    "bivariatePinkGreen",
    "bivariateRedCyan",
    "bivariateGreenRed",
    "bivariateBlueYellow"
  )

  palettes <- listPalettes()
  expect_s3_class(palettes, "data.frame")
  expect_true(all(c("name", "type", "n", "colours", "is_default") %in% colnames(palettes)))
  expect_true(all(c("categorical", "continuous", "bivariate") %in% palettes$type))
  expect_false("discrete" %in% palettes$type)
  expect_false("default" %in% palettes$name)
  expect_true(all(expected_bivariate %in% palettes$name))
  expect_s3_class(palettes$colours, "AsIs")

  categorical <- listPalettes(type = "categorical")
  expect_true(all(categorical$type == "categorical"))
  expect_true("soft" %in% categorical$name)
  expect_true("combinationBright" %in% categorical$name)
  combination_row <- which(categorical$name == "combinationBright")
  expect_equal(categorical$n[combination_row], 32L)
  expect_equal(length(categorical$colours[[combination_row]]), 32L)
  expect_true(categorical$is_default[categorical$name == "scOverlay"])
  expect_true(all(!categorical$is_default[categorical$name != "scOverlay"]))

  continuous <- listPalettes(type = "continuous")
  expect_true(all(continuous$type == "continuous"))
  expect_true("gray_red" %in% continuous$name)
  gray_red_row <- which(continuous$name == "gray_red")
  viridis_row <- which(continuous$name == "viridis")
  expect_equal(continuous$n[gray_red_row], 5L)
  expect_equal(length(continuous$colours[[gray_red_row]]), 5L)
  expect_true("viridis" %in% continuous$name)
  expect_true(is.na(continuous$n[viridis_row]))
  expect_true(is.na(continuous$colours[[viridis_row]]))
  expect_true(continuous$is_default[continuous$name == "gray_red"])
  expect_true(all(!continuous$is_default[continuous$name != "gray_red"]))

  bivariate <- listPalettes(type = "bivariate")
  expect_true(all(bivariate$type == "bivariate"))
  expect_true(all(expected_bivariate %in% bivariate$name))
  expect_true(bivariate$is_default[bivariate$name == "bivariateBright"])
  expect_true(all(!bivariate$is_default[bivariate$name != "bivariateBright"]))
  expect_equal(bivariate$n, rep(4L, nrow(bivariate)))
  bivariate_bright_row <- which(bivariate$name == "bivariateBright")
  expect_named(
    bivariate$colours[[bivariate_bright_row]],
    c("low_xy", "high_x", "high_y", "high_xy")
  )

  expect_error(listPalettes(type = "discrete"), "`type` must be one of")
})

test_that("plotPalettes can draw palette previews", {
  expect_error(listPalettes(plot = TRUE), "unused argument")

  p_bivariate <- plotPalettes(type = "bivariate")
  expect_true(inherits(p_bivariate, "ggplot") || inherits(p_bivariate, "patchwork"))

  p_all <- plotPalettes(type = "all")
  expect_true(inherits(p_all, "ggplot") || inherits(p_all, "patchwork"))

  expect_error(plotPalettes(type = "discrete"), "`type` must be one of")

  metadata <- listPalettes()
  expect_s3_class(metadata, "data.frame")
  expect_named(metadata, c("name", "type", "n", "colours", "is_default"))
})

test_that("plotOverlay works with default and named continuous foreground palettes", {
  sce <- get_example_sce()

  p_default <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    fg_type = "continuous"
  )
  expect_s3_class(p_default, "ggplot")

  p_gray_red <- plotOverlay(
    sce,
    reduced_dim = "TSNE",
    background = "cluster",
    foreground = "SOX10",
    fg_type = "continuous",
    fg_palette = "gray_red"
  )
  expect_s3_class(p_gray_red, "ggplot")
})
