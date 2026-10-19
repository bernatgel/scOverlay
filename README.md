<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/bernatgel/scOverlay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bernatgel/scOverlay/actions/workflows/R-CMD-check.yaml)
[![Bioconductor status](https://bioconductor.org/shields/build/release/bioc/scOverlay.svg)](https://bioconductor.org/packages/scOverlay) -->
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/license/artistic-2-0)
<!-- badges: end -->

# scOverlay: Multilayer plots for single-cell data

<p align="center">
  <img width="6000" height="1500" alt="scOverlay example plot" src="https://github.com/user-attachments/assets/8403e754-a709-4e9d-80ad-ced2aa681fbb"/>
</p>


## Description

`scOverlay` creates layered visualizations of single-cell data stored in
`SingleCellExperiment` objects. It is designed to show one layer of information
as context, such as clusters, samples or gene expression, while overlaying a
second layer of information on top of it.

The package is especially useful to highlight subsets of cells, gene expression,
metadata annotations or quality-control scores while preserving the global
structure of a t-SNE, UMAP or other reduced-dimension representation.

## Installation

`scOverlay` is a Bioconductor package (in process). Once available from Bioconductor, it
can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("scOverlay")
```

The development version can be installed from GitHub with:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("bernatgel/scOverlay")
```

## Getting started

```r
library(scOverlay)

data("sce_overlay_example")

plotOverlay(
    sce = sce_overlay_example,
    foreground = "SOX10",
    background = "cluster",
    reduced_dim = "TSNE"
)
```

`background` and `foreground` can be a gene, a cell annotation
from `colData(sce)`, `"solid"` or `"none"`.

```r
plotOverlay(
    sce = sce_overlay_example,
    foreground = "S100B",
    background = "cluster",
    fg_subset_cells = quote(cluster == "C1"),
    fg_palette = "magma",
    reduced_dim = "TSNE"
)
```

Grouped plots can be created with `plotOverlayPerGroup()`:

```r
plotOverlayPerGroup(
    sce = sce_overlay_example,
    group_col = "sample_type",
    foreground = "SOX10",
    background = "cluster",
    reduced_dim = "TSNE",
    fg_limits = "shared",
    shared_legend = TRUE
)
```

## Documentation

The package vignette provides a more complete introduction to the layered
plotting model, palettes, foreground subsetting, grouped plots, rasterization
and saving figures.

- Bioconductor landing page: <https://bioconductor.org/packages/scOverlay>
- Vignette: <https://bioconductor.org/packages/release/bioc/vignettes/scOverlay/inst/doc/scOverlay.html>
- Development repository: <https://github.com/bernatgel/scOverlay> 

Bioconductor links will become active once the package is available from Bioconductor.

## Main features

- Dual-layer embedding plots with independent background and foreground colour scales.
- Foreground filtering while retaining the full embedding as context.
- Grouped and nested plots with compact or grid layouts and optional shared scales and legends.
- Bivariate colouring of two numeric foreground features.
- Combination colouring of two to five numeric foreground features.
- Palette listing, inspection and plotting utilities.
- Optional rasterization of dense background and foreground point layers.
- Standard `ggplot2` output objects that can be further customized.

## Citation

If you use `scOverlay` in published work, please cite it with:

```r
citation("scOverlay")
```


## License

`scOverlay` is distributed under the Artistic-2.0 license.

## Bug reports and feature requests

Please use the GitHub issue tracker to report bugs or request features:

<https://github.com/bernatgel/scOverlay/issues>
