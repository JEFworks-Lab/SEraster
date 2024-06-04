# Spatial Experiments raster (SEraster)

[![R-CMD-check](https://github.com/JEFworks-Lab/SEraster/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/JEFworks-Lab/SEraster/actions/workflows/check-standard.yaml)

`SEraster` is a rasterization preprocessing framework that aggregates cellular information into spatial pixels to reduce resource requirements for spatial omics data analysis. This is the `SEraster` R documentation website. Questions, suggestions, or problems should be submitted as [GitHub issues](https://github.com/JEFworks-Lab/SEraster/issues).

<p>

<img src="https://github.com/JEFworks/SEraster/blob/main/images/seraster_logo_hex.png?raw=true" align="center" height="300" style="float: center; height:300px;"/>

</p>

## Overview

`SEraster` reduces the number of spatial points in spatial omics datasets for downstream analysis through a process of rasterization where single cells' gene expression or cell-type labels are aggregated into equally sized pixels based on a user-defined `resolution`. Here, we refer to a particular `resolution` of rasterization by the side length of the pixel such that finer `resolution` indicates smaller pixel size and coarser `resolution` indicates larger pixel size.

<p align="center">

<img src="https://github.com/JEFworks-Lab/SEraster/blob/main/images/overview.png?raw=true" height="600"/>

</p>

## Installation

To install `SEraster`, we currently recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/SEraster')
```

## Tutorials

Introduction:

-   [Formatting a SpatialExperiment Object for SEraster](https://jef.works/SEraster/articles/formatting-SpatialExperiment-for-SEraster.html)
-   [Getting Started With SEraster](https://jef.works/SEraster/articles/getting-started-with-SEraster.html)

## Citation

Our preprint describing `SEraster` is available on *bioRxiv*:

[Aihara G. et al. (2024), "SEraster: a rasterization preprocessing framework for scalable spatial omics data analysis", *bioRxiv*](https://doi.org/10.1101/2024.02.01.578436)