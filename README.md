# Spatial Experiments raster (SEraster)

`SEraster` is a pre-processing tool to enable scalable and accurate analysis of large-scale spatial omics datasets with existing tools.

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
``` r

```

### Downstream Analysis
``` r

```

#### Spatial variable gene (SVG) analysis
``` r

```

#### Cell-type cooccurrence analysis
``` r

```

## Citation
