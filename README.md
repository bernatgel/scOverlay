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

`scOverlay` is currently under review for inclusion in Bioconductor. Once
available from Bioconductor, it can be installed with:

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

data("sce_mpnst_example")

plotOverlay(
    sce = sce_mpnst_example,
    foreground = "SOX10",
    background = "cluster",
    reduced_dim = "TSNE",
    bg_legend = FALSE
)

```
<img width="2100" height="1800" alt="SOX10 expression over cell clusters" src="https://github.com/user-attachments/assets/e75a38f9-d818-4dd5-9338-0a9d0bd2bbf4" />


`background` and `foreground` can be a gene, a cell annotation
from the SingleCellExperiment object (`colData(sce)`)
or a solid colour.

We can filter the cells in the foreground using any criterion,
including cell annotations in `colData(sce)`.


```r
plotOverlay(
    sce = sce_mpnst_example,
    foreground = "CDH19",
    background = "cluster",
    fg_subset_cells = quote(celltype.main == "Schwann cell"),
    fg_palette = "viridis",
    reduced_dim = "TSNE",
    bg_legend = FALSE
)

```
<img width="2100" height="1800" alt="CDH19 expression in Schwann cells" src="https://github.com/user-attachments/assets/efedb7d2-ecb0-482c-bae1-99b8887f5015" />


`plotOverlayPerGroup()` groups plots by a metadata variable:

```r
plotOverlayPerGroup(
    sce = sce_mpnst_example,
    group_col = "sample_type",
    foreground = "CDH19",
    background = "celltype.main",
    reduced_dim = "TSNE",
    fg_limits = "shared",
    shared_legend = TRUE
)
```
<img width="4800" height="1200" alt="CDH19 overlays grouped by sample type" src="https://github.com/user-attachments/assets/8bceaaf0-a8c5-4b5f-a711-391f2c8c797d" />


scOverlay supports bivariate colouring for two numeric variables and combination
colouring for two to five numeric variables. The example below uses combination
colouring for three genes:


```r
plotOverlay(
    sce = sce_mpnst_example,
    foreground = c("PRRX1", "EGFR", "CDH19"),
    fg_type = "combination",
    background = "celltype.main",
    reduced_dim = "TSNE"
)
```
<img width="2400" height="1800" alt="Combination colouring of three genes" src="https://github.com/user-attachments/assets/e13c1664-09fb-48ef-bcc2-ec9997b7cf48" />


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
