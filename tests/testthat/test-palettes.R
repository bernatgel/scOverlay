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
  palettes <- listPalettes()
  expect_s3_class(palettes, "data.frame")
  expect_true(all(c("name", "type", "n_colours", "is_default") %in% colnames(palettes)))
  expect_true(all(c("categorical", "continuous") %in% palettes$type))
  expect_false("discrete" %in% palettes$type)
  expect_false("default" %in% palettes$name)

  categorical <- listPalettes(type = "categorical")
  expect_true(all(categorical$type == "categorical"))
  expect_true("soft" %in% categorical$name)
  expect_true("combinationBright" %in% categorical$name)
  expect_equal(categorical$n_colours[categorical$name == "combinationBright"], 32L)
  expect_true(categorical$is_default[categorical$name == "scOverlay"])
  expect_true(all(!categorical$is_default[categorical$name != "scOverlay"]))

  continuous <- listPalettes(type = "continuous")
  expect_true(all(continuous$type == "continuous"))
  expect_true("gray_red" %in% continuous$name)
  expect_true(continuous$is_default[continuous$name == "gray_red"])
  expect_true(all(!continuous$is_default[continuous$name != "gray_red"]))

  expect_error(listPalettes(type = "discrete"), "'arg' should be one of")
})

test_that("listPalettes can draw a ggplot swatch preview", {
  p <- listPalettes(plot = TRUE, n = 16)
  expect_s3_class(p, "ggplot")

  p_names <- listPalettes(plot = TRUE, n = 16, show_names = TRUE)
  expect_s3_class(p_names, "ggplot")

  p_no_names <- listPalettes(plot = TRUE, n = 16, show_names = FALSE)
  expect_s3_class(p_no_names, "ggplot")

  metadata <- listPalettes(plot = FALSE)
  expect_s3_class(metadata, "data.frame")
  expect_named(metadata, c("name", "type", "n_colours", "is_default"))
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
