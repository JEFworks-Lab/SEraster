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

To install `SEraster` using Bioconductor, start R (version "4.4.0") and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SEraster")
```

The latest development version can also be installed from [GitHub](https://github.com/JEFworks-Lab/SEraster) using `remotes`:

```r
require(remotes)
remotes::install_github('JEFworks-Lab/SEraster')
```

In addition, `SEraster` is also compatible with `SeuratObject` through `SeuratWrappers`. `SeuratWrappers` implementation can be installed using `remotes`:

```r
require(remotes)
remotes::install_github('satijalab/seurat-wrappers@SEraster')
```

Documentation and tutorial for the `SeuratWrappers` implementation can be found in the `SEraster` branch of the [`SeuratWrappers` GitHub repository](https://github.com/satijalab/seurat-wrappers/tree/SEraster).

## Tutorials

Introduction:

-   [Formatting a SpatialExperiment Object for SEraster](https://jef.works/SEraster/articles/formatting-SpatialExperiment-for-SEraster.html)
-   [Getting Started With SEraster](https://jef.works/SEraster/articles/getting-started-with-SEraster.html)
-   [SEraster for Spatial Variable Genes Analysis](https://jef.works/SEraster/articles/SEraster-for-SVG-analysis.html)
-   [Characterizing mPOA cell-type heterogeneity with spatial bootstrapping](https://jef.works/SEraster/articles/characterizing-mPOA-cell-type-heterogeneity.html)

## Citation

Our manuscript describing `SEraster` is available on *Bioinformatics*:

[Gohta Aihara, Kalen Clifton, Mayling Chen, Zhuoyan Li, Lyla Atta, Brendan F Miller, Rahul Satija, John W Hickey, Jean Fan, SEraster: a rasterization preprocessing framework for scalable spatial omics data analysis, Bioinformatics, Volume 40, Issue 7, July 2024, btae412, https://doi.org/10.1093/bioinformatics/btae412](https://academic.oup.com/bioinformatics/article/40/7/btae412/7696710)