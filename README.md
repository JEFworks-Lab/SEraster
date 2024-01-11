# Spatial Experiments raster (SEraster)

`SEraster` is a pre-processing tool to enable scalable and accurate analysis of large-scale spatial omics datasets with existing tools.

<p align="center">
  <img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/seraster_logo_hex.png?raw=true" height="200"/>
</p>

## Overview

## Installation

To install `SEraster`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/SEraster')
```

## Tutorials

## Input data format

In the examples below, we assume the input data is provided as a `SpatialExperiment` Bioconductor object. Please refer to the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) package and the `merfish_mousePOA` dataset in the package to see how you would format your data into a `SpatialExperiment` object.

## Example

A short example workflow is shown below.

### Load packages
``` r
library(SpatialExperiment)
library(SEraster)
```

### Load example dataset
``` r
data("merfish_mousePOA")
dim(merfish_mousePOA)
```

### Getting started

#### Rasterize gene expression
``` r
rastGexp <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 50)

# plot total rasterized gene expression
SEraster::plotRaster(rastGexp, name = "Total rasterized gene expression")
```
<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_gexp_total.png?raw=true" height="500"/>
</p>

``` r
# plot specific gene
SEraster::plotRaster(rastGexp, feature_name = "Esr1", name = "Esr1")
```

#### Rasterize cell type labels
``` r
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50)

# plot
SEraster::plotRaster(rastCt, name = "cell counts", option = "inferno")
SEraster::plotRaster(rastCt, feature_name = "Inhibitory", name = "Inhibitory neuron counts", option = "inferno")
```

### Downstream Analysis

We will highlight spatial variable gene (SVG) and cell-type cooccurrence analyses as examples of downstream analysis that benefit from SEraster.

#### Spatial variable gene (SVG) analysis

Here, we use a previously developed tool called `nnSVG`. Please refer to [nnSVG](https://bioconductor.org/packages/nnSVG) for more details about the package. We can directly input rasterized gene expression `SpatialExperiment` object from `SEraster` into `nnSVG`.

``` r
# run nnSVG
set.seed(0)
rastGexp <- nnSVG(rastGexp, assay_name = "pixelval")
```

``` r
# number of significant SVGs
table(rowData(rastGexp)$padj <= 0.05)
```

``` r
# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp)$rank == 1)
top_svg_name <- rowData(rastGexp)$gene_name[top_svg]
SEraster::plotRaster(rastGexp, feature_name = top_svg_name, name = top_svg_name)
```

We can also perform cell-type specific SVG analysis by subsetting the dataset prior to applying SEraster.

``` r
# subset data
ct_interest <- "Ependymal"
spe_epend <- merfish_mousePOA[,merfish_mousePOA$celltype == ct_interest]

# run SEraster
rastGexp_epend <- SEraster::rasterizeGeneExpression(spe_epend, assay_name="volnorm", resolution = 50)

# run nnSVG
set.seed(0)
rastGexp_epend <- nnSVG(rastGexp_epend, assay_name = "pixelval")
```

``` r
# number of significant SVGs
table(rowData(rastGexp_epend)$padj <= 0.05)
```

``` r
# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp_epend)$rank == 1)
top_svg_name <- rowData(rastGexp_epend)$gene_name[top_svg]
SEraster::plotRaster(rastGexp_epend, feature_name = top_svg_name, name = top_svg_name)
```

#### Cell-type cooccurrence analysis

Rasterized cell-type labels can be used to analyze pair-wise cell-type cooccurrence. To do so, we binarize the rasterized cell-type labels using a relative enrichment metric and a previously developed tool called `CooccurrenceAffinity`. Please refer to our paper for more details about the methodology and [CooccurrenceAffinity](https://CRAN.R-project.org/package=CooccurrenceAffinity) for more details about the package.

``` r
# extract cell-type labels
ct_labels <- as.factor(colData(merfish_mousePOA)$celltype)

# compute relative enrichment (RE) metric
mat <- assay(rastCt, "pixelval")
mat_re <- do.call(rbind, lapply(rownames(rastCt), function(ct_label) {
    mat[ct_label,] / (sum(mat[ct_label,]) / sum(mat) * colSums(mat))
}))
rownames(mat_re) <- rownames(mat)

# binarize
mat_bin <- ifelse(mat_re >= 1, 1, 0)

# add RE and binarized layers to SpatialExperiment object
assays(rastCt) <- list(pixelval = assay(rastCt, "pixelval"), re = mat_re, bin = mat_bin)
```

``` r
ct_interest <- "Ependymal"

# plot pixel value for a cell-type of interest
plotRaster(rastCt, assay_name = "pixelval", feature_name = ct_interest, name = "cell-type counts", option = "inferno")
```

``` r
# plot RE value for a cell-type of interest
plotRaster(rastCt, assay_name = "re", feature_name = ct_interest, name = "RE", option = "inferno")
```

``` r
# plot binarized value for a cell-type of interest
plotRaster(rastCt, assay_name = "bin", feature_name = ct_interest, factor_levels = c(0,1), name = "binarized", option = "inferno")
```

``` r
# run CooccurrenceAffinity
ct_coocc <- CooccurrenceAffinity::affinity(data = mat_bin, row.or.col = "row", squarematrix = c("all"))

# plot maximum likelihood estimates of affinity metric (alpha MLE)
CooccurrenceAffinity::plotgg(data = ct_coocc, variable = "alpha_mle", legendlimit = "datarange")
```

## Citation
