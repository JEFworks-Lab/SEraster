# Spatial Experiments raster (SEraster)

`SEraster` is a rasterization preprocessing framework that aggregates cellular information into spatial pixels to reduce resource requirements for spatial omics data analysis.

<p align="center">
  <img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/seraster_logo_hex.png?raw=true" height="200"/>
</p>

## Overview

`SEraster` reduces the number of spatial points in spatial omics datasets for downstream analysis through a process of rasterization where single cells’ gene expression or cell-type labels are aggregated into equally sized pixels based on a user-defined `resolution`. Here, we refer to a particular `resolution` of rasterization by the side length of the pixel such that finer `resolution` indicates smaller pixel size and coarser `resolution` indicates larger pixel size.

<p align="center">
  <img src="https://github.com/JEFworks-Lab/SEraster/blob/main/docs/images/overview.png?raw=true" height="600"/>
</p>

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

# check the dimension of the genes-by-cells matrix at single-cell resolution
dim(merfish_mousePOA)
```

``` r
[1]  155 6509
```

``` r
# check the number of cell-types
length(unique(colData(merfish_mousePOA)$celltype))
```

``` r
[1]  16
```

This MERFISH mouse preoptic area dataset contains 6,509 cells and 16 cell-types.

``` r
# plot at single-cell resolution
df <- data.frame(spatialCoords(merfish_mousePOA), celltype = colData(merfish_mousePOA)$celltype)
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 1) +
  labs(x = "x (μm)",
       y = "y (μm)",
       col = "Cell-types") +
  theme_bw() +
  theme(panel.grid = element_blank())
```

<p align="center">
  <img src="https://github.com/JEFworks-Lab/SEraster/blob/main/docs/images/singlecell_celltypes.png?raw=true" height="550"/>
</p>

### Getting started

#### Rasterize gene expression
``` r
rastGexp <- SEraster::rasterizeGeneExpression(merfish_mousePOA, assay_name="volnorm", resolution = 50)

# check the dimension of the genes-by-cells matrix after rasterizing gene expression
dim(rastGexp)
```
``` r
[1]  155 1301
```
As you can see, SEraster aggregated 6,509 single cells into 1,301 pixels.

``` r
# plot total rasterized gene expression
SEraster::plotRaster(rastGexp, name = "Total rasterized gene expression")
```
<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_gexp_total.png?raw=true" height="400"/>
</p>

``` r
# plot specific gene
SEraster::plotRaster(rastGexp, feature_name = "Esr1", name = "Esr1")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_gexp_esr1.png?raw=true" height="400"/>
</p>

#### Rasterize cell-type labels
``` r
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50)

# check the dimension of the cell-types-by-cells matrix after rasterizing cell-type labels
dim(rastGexp)
```
``` r
[1]  16 1301
```
``` r
# plot total cell counts
SEraster::plotRaster(rastCt, name = "cell counts", option = "inferno")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_ct_total.png?raw=true" height="400"/>
</p>

``` r
# plot specific cell-type
SEraster::plotRaster(rastCt, feature_name = "Inhibitory", name = "Inhibitory neuron counts", option = "inferno")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_ct_inhibitory.png?raw=true" height="400"/>
</p>

### Sample Downstream Analysis

`SEraster` returns rasteriezd gene-expression and cell-type information as `SpatialExperiment` objects that can be integrated with other existing downstream analysis tools. We demonstrate below spatial variable gene (SVG) and cell-type cooccurrence analyses as examples of such potential downstream analysis.

#### Spatial variable gene (SVG) analysis

Here, we use a previously developed tool called `nnSVG`. Please refer to [nnSVG](https://bioconductor.org/packages/nnSVG) for more details about the package. We can directly input rasterized gene expression `SpatialExperiment` object from `SEraster` into `nnSVG`.

``` r
library(nnSVG)
```

``` r
# run nnSVG
set.seed(0)
rastGexp <- nnSVG(rastGexp, assay_name = "pixelval")
```

``` r
# number of significant SVGs
table(rowData(rastGexp)$padj <= 0.05)
```

```r
## 
## FALSE  TRUE
##    17   138
```

``` r
# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp)$rank == 1)
top_svg_name <- rownames(rowData(rastGexp))[top_svg]
SEraster::plotRaster(rastGexp, feature_name = top_svg_name, name = top_svg_name)
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_gexp_top_svg.png?raw=true" height="400"/>
</p>

We can also perform cell-type specific SVG analysis by subsetting the dataset prior to applying SEraster.

``` r
# subset data
ct_interest <- "Excitatory"
spe_sub <- merfish_mousePOA[,merfish_mousePOA$celltype == ct_interest]

# run SEraster
rastGexp_sub <- SEraster::rasterizeGeneExpression(spe_sub, assay_name="volnorm", resolution = 50)

# run nnSVG
set.seed(0)
rastGexp_sub <- nnSVG(rastGexp_sub, assay_name = "pixelval")
```

``` r
# number of significant SVGs
table(rowData(rastGexp_sub)$padj <= 0.05)
```

``` r
## 
## FALSE  TRUE 
##    45   110
```

``` r
# plot rasterized gene expression of top-ranked SVG
top_svg <- which(rowData(rastGexp_sub)$rank == 1)
top_svg_name <- rownames(rowData(rastGexp_sub))[top_svg]
SEraster::plotRaster(rastGexp_sub, feature_name = top_svg_name, name = top_svg_name)
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_gexp_sub_top_svg.png?raw=true" height="400"/>
</p>

#### Cell-type cooccurrence analysis

Rasterized cell-type labels can be used to analyze pair-wise cell-type cooccurrence. To do so, we binarize the rasterized cell-type labels using a relative enrichment metric and a previously developed tool called `CooccurrenceAffinity`. Please refer to our paper for more details about the methodology and [CooccurrenceAffinity](https://CRAN.R-project.org/package=CooccurrenceAffinity) for more details about the package.

``` r
library(CooccurrenceAffinity)
```

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

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_ct_sub_total.png?raw=true" height="400"/>
</p>

``` r
# plot RE value for a cell-type of interest
plotRaster(rastCt, assay_name = "re", feature_name = ct_interest, name = "RE", option = "inferno")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_ct_sub_re.png?raw=true" height="400"/>
</p>

``` r
# plot binarized value for a cell-type of interest
plotRaster(rastCt, assay_name = "bin", feature_name = ct_interest, factor_levels = c(0,1), name = "binarized", option = "inferno")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/rasterized_ct_sub_bin.png?raw=true" height="400"/>
</p>

``` r
# run CooccurrenceAffinity
ct_coocc <- CooccurrenceAffinity::affinity(data = mat_bin, row.or.col = "row", squarematrix = c("all"))

# plot maximum likelihood estimates of affinity metric (alpha MLE)
CooccurrenceAffinity::plotgg(data = ct_coocc, variable = "alpha_mle", legendlimit = "datarange")
```

<p align="center">
<img src="https://github.com/JEFworks/SEraster/blob/main/docs/images/coocc_heatmap.png?raw=true" height="500"/>
</p>

## Citation

Our preprint describing `SEraster` is available on bioRxiv:

[Aihara G. et al. (2024), "SEraster: a rasterization preprocessing framework for scalable spatial omics data analysis", bioRxiv](https://doi.org/10.1101/2024.02.01.578436)
