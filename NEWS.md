# scOverlay 0.99.5

## Initial Bioconductor submission

* Introduced `plotOverlay()` for displaying a foreground feature over a
  contextual background in a two-dimensional `SingleCellExperiment` embedding.
  Background and foreground layers have independent colour scales and can
  represent gene expression, cell metadata, or a solid colour.
* Added foreground-only cell selection, preserving the full embedding as
  background context. Supports configurable assays, point sizes, transparency,
  colour limits, legends, and drawing order.
* Added `plotGeneOverlay()` for plotting one or more gene-expression overlays.
* Added `plotOverlayPerGroup()` and `plotOverlayPerNestedGroup()` for comparing
  subsets across grouped panels, with optional shared foreground scales and
  legends. Nested layouts support automatic selection between compact and
  aligned grids, with informative panel titles.
* Added bivariate foreground colouring for two numeric features, with
  `plotBivariateLegend()` for the corresponding two-dimensional legend.
* Added combination colouring for the thresholded positive/negative patterns
  of two to five numeric features.
* Added palette utilities `getPalette()`, `listPalettes()`, and `plotPalettes()`
  for selecting, inspecting, and previewing categorical, continuous, and
  bivariate colour palettes.
* Added optional point rasterization for dense embeddings while retaining
  `ggplot2`-compatible plot output.
