#' rasterizeSparseMatrix
#' 
#' @description Function to rasterize a given input sparse matrix based on a given position matrix.
#' 
#' @description This function assumes that inputs are provided as a \code{dgCmatrix} object 
#' and \code{matrix} for data and position, respectively.
#' 
#' @param data \code{dgCmatrix}: Feature x observation matrix represented as a 
#' \code{dgCmatrix} object. Features can be genes or cell types. In the case of 
#' features being cell types, this matrix is assumed to be a sparse model matrix 
#' with rows as cell types and columns as cell IDs.
#' 
#' @param pos \code{matrix}: Spatial x,y coordinates of observations stored as a 
#' matrix array. Further, x,y coordinates are assumed to be stored in column 1 
#' and 2 of \code{spatialCoords}.
#' 
#' @param resolution \code{integer}: Resolution or side length of each pixel. 
#' The unit of this parameter is assumed to be the same as the unit of spatial 
#' coordinates of the input data.
#' 
#' @param fun \code{character}: If "mean", pixel value for each pixel 
#' would be the average of gene expression for all cells within the pixel. If 
#' "sum", pixel value for each pixel would be the sum of gene expression for all 
#' cells within the pixel.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization. Default = 1. 
#' Inputting this argument when the \code{BPPARAM} argument is missing would set parallel 
#' exeuction back-end to be \code{BiocParallel::MulticoreParam(workers = n_threads)}. 
#' We recommend setting this argument to be the number of cores available 
#' (\code{parallel::detectCores(logical = FALSE)}). If \code{BPPARAM} argument is 
#' not missing, the \code{BPPARAM} argument would override \code{n_threads} argument.
#' 
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for parallelization. 
#' This argument is provided for advanced users of \code{BiocParallel} for further 
#' flexibility for setting up parallel-execution back-end. Default is NULL. If 
#' provided, this is assumed to be an instance of \code{BiocParallelParam}.
#' 
#' @return The output is returned as a \code{list} containing rasterized feature 
#' x observation matrix as \code{dgCmatrix}, spatial x,y coordinates of pixel 
#' centroids as \code{matrix}, and \code{data.frame} containing cell IDs of cells 
#' that were aggregated in each pixel.
#' 
#' @importFrom raster extent res rasterToPoints cellFromXY
#' @importFrom methods new
#' @importFrom BiocParallel MulticoreParam bpstart bplapply
#' @importFrom Matrix rowMeans rowSums
#' 
#' @export
#' 
rasterizeSparseMatrix <- function(data, pos, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL) {
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_threads)
  }
  BiocParallel::bpstart(BPPARAM)
  
  ## create RasterLayer (simplified way, using actual operations in raster() function)
  ext <- raster::extent(min(pos[,1])-resolution/2, max(pos[,1])+resolution/2, min(pos[,2])-resolution/2, max(pos[,2])+resolution/2)
  r <- methods::new('RasterLayer', extent=ext)
  raster::res(r) <- resolution
  
  ## convert RasterLayer to SpatialPoints (sp package)
  ## might be good to use sf package since sp package is retiring soon
  pts <- raster::rasterToPoints(r, spatial = TRUE)
  
  ## get x,y coordinates for each pixel
  pos_pixel <- pts@coords
  rownames(pos_pixel) <- paste0("pixel", seq_along(pos_pixel[,1]))
  
  ## assign pixel ID to each single cell xy coordinates
  pixel_ids <- raster::cellFromXY(r, pos)
  names(pixel_ids) <- rownames(pos)
  
  ## store aggregated subsetted sparse matrix data for each pixel into a CsparseMatrix, store cell IDs and # of cells for each pixel into a data frame
  out <- unlist(BiocParallel::bplapply(seq_along(pos_pixel[,1]), function(id){
    ## get cell IDs for a particular pixel
    cell_ids <- names(pixel_ids[pixel_ids == id])
    ## subset feature observation matrix
    spmat <- data[,cell_ids, drop = FALSE]
    ## aggregate cell counts to create pixel value
    if (fun == "mean") {
      pixel_val <- Matrix::rowMeans(spmat, sparseResult = TRUE) ## returns numeric if sparseResult = FALSE, dsparseVector if sparseResult = TRUE
    } else if (fun == "sum") {
      pixel_val <- Matrix::rowSums(spmat, sparseResult = TRUE) ## returns numeric if sparseResult = FALSE, dsparseVector if sparseResult = TRUE
    }
    
    meta_rast <- data.frame(num_cell = length(cell_ids))
    meta_rast$cellID_list <- list(cell_ids)
    
    return(list(as(pixel_val, "CsparseMatrix"), meta_rast))
  }), recursive = FALSE)
  
  ## extract rasterized sparse matrix
  data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
  
  ## extract rasterized data frame
  meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
  
  ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
  rownames(data_rast) <- rownames(data)
  colnames(data_rast) <- paste0("pixel", seq_along(pos_pixel[,1]))
  rownames(meta_rast) <- paste0("pixel", seq_along(pos_pixel[,1]))
  
  ## convert NaN to NA
  data_rast[is.nan(data_rast)] <- NA
  
  ## output
  output <- list("data_rast" = data_rast, "pos_rast" = pos_pixel, "meta_rast" = meta_rast)
}

#' rasterizeGeneExpression
#' 
#' @description Function to rasterize feature x observation matrix in spatially-resolved 
#' transcriptomics data represented as SpatialExperiment class.
#'  
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} object.
#' 
#' @param input \code{SpatialExperiment}: Input data represented as a \code{SpatialExperiment} object. 
#' It is assumed to have an \code{assay} slot containing feature (genes) x observation (cells)
#' matrix as \code{dgCmatrix} and a \code{spatialCoords} slot containing spatial x,y 
#' coordinates of observations as matrix array. Further, x,y coordinates are assumed 
#' to be stored in coloumn 1 and 2 of \code{spatialCoords}.
#' 
#' @param assay_name \code{character}: Name of the assay slot of the input that 
#' you want to apply rasterization. If no argument is given, the first assay of the input 
#' would be rasterized. This argument is useful when you have both raw and normalized 
#' assays stored in the input, and you want to apply rasterization to the normalized assay.
#' 
#' @param resolution \code{integer}: Resolution or side length of each pixel. 
#' The unit of this parameter is assumed to be the same as the unit of spatial 
#' coordinates of the input data.
#' 
#' @param fun \code{character}: If "mean", pixel value for each pixel 
#' would be the average of gene expression for all cells within the pixel. If 
#' "sum", pixel value for each pixel would be the sum of gene expression for all 
#' cells within the pixel.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization. Default = 1. 
#' Inputting this argument when the \code{BPPARAM} argument is missing would set parallel 
#' exeuction back-end to be \code{BiocParallel::MulticoreParam(workers = n_threads)}. 
#' We recommend setting this argument to be the number of cores available 
#' (\code{parallel::detectCores(logical = FALSE)}). If \code{BPPARAM} argument is 
#' not missing, the \code{BPPARAM} argument would override \code{n_threads} argument.
#' 
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for parallelization. 
#' This argument is provided for advanced users of \code{BiocParallel} for further 
#' flexibility for setting up parallel-execution back-end. Default is NULL. If 
#' provided, this is assumed to be an instance of \code{BiocParallelParam}.
#' 
#' @param na.rm \code{logical}: If TRUE, NA values in the output 
#' rasterized feature (genes) x observation (pixels) matrix (dgCmatrix) and 
#' spatial coordinates matrix (matrix array) are removed. This corresponds to 
#' removing pixels with no cells.
#' 
#' @return The output is returned as a new \code{SpatialExperiment} object with 
#' \code{assay} slot containing the feature (genes) x observations (pixels) matrix 
#' (dgCmatrix), \code{spatialCoords} slot containing spatial x,y coordinates of 
#' pixel centroids, and \code{colData} slot containing cell IDs of cells that 
#' were aggregated in each pixel.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment assay
#' 
#' @export
#' 
rasterizeGeneExpression <- function(input, assay_name = NULL, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL, na.rm = FALSE) {
  ## rasterize
  if (is.null(assay_name)) {
    out <- rasterizeSparseMatrix(assay(input), spatialCoords(input), resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
  } else {
    stopifnot(is.character(assay_name))
    out <- rasterizeSparseMatrix(assay(input, assay_name), spatialCoords(input), resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
  }
  data_rast <- out$data_rast
  pos_rast <- out$pos_rast
  meta_rast <- out$meta_rast
  
  ## remove NA based on the na.rm argument
  if (na.rm) {
    na_cols <- colSums(is.na(data_rast)) != 0
    data_rast <- data_rast[,!na_cols]
    pos_rast <- pos_rast[!na_cols,]
    meta_rast <- meta_rast[!na_cols,]
  }
  
  ## construct a new SpatialExperiment object as output
  output <- SpatialExperiment::SpatialExperiment(
    assays = list(pixelval = data_rast),
    spatialCoords = pos_rast,
    colData = meta_rast
  )
  
  return(output)
}

#' rasterizeCellType
#' 
#' @description Function to rasterize cell type labels in spatially-resolved 
#' transcriptomics data represented as SpatialExperiment class.
#'  
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} object.
#' 
#' @param input \code{SpatialExperiment}: Input data represented as a \code{SpatialExperiment} object. 
#' It is assumed to have a \code{colData} slot containing cell type 
#' labels for observations as a data frame column and a \code{spatialCoords} slot 
#' containing spatial x,y coordinates of observations as matrix array. Further, 
#' x,y coordinates are assumed to be stored in coloumn 1 and 2 of \code{spatialCoords}.
#' 
#' @param col_name \code{character}: Column name of the \code{colData} object 
#' containing cell type labels for observations.
#' 
#' @param resolution \code{integer}: Resolution or side length of each pixel. 
#' The unit of this parameter is assumed to be the same as the unit of spatial 
#' coordinates of the input data.
#' 
#' @param fun \code{character}: If "mean", pixel value for each pixel 
#' would be the average of gene expression for all cells within the pixel. If 
#' "sum", pixel value for each pixel would be the sum of gene expression for all 
#' cells within the pixel.
#' 
#' @param n_threads \code{integer}: Number of threads for parallelization. Default = 1. 
#' Inputting this argument when the \code{BPPARAM} argument is missing would set parallel 
#' exeuction back-end to be \code{BiocParallel::MulticoreParam(workers = n_threads)}. 
#' We recommend setting this argument to be the number of cores available 
#' (\code{parallel::detectCores(logical = FALSE)}). If \code{BPPARAM} argument is 
#' not missing, the \code{BPPARAM} argument would override \code{n_threads} argument.
#' 
#' @param BPPARAM \code{BiocParallelParam}: Optional additional argument for parallelization. 
#' This argument is provided for advanced users of \code{BiocParallel} for further 
#' flexibility for setting up parallel-execution back-end. Default is NULL. If 
#' provided, this is assumed to be an instance of \code{BiocParallelParam}.
#' 
#' @param na.rm \code{logical}: If TRUE, NA values in the output 
#' rasterized feature (genes) x observation (pixels) matrix (dgCmatrix) and 
#' spatial coordinates matrix (matrix array) are removed. This corresponds to 
#' removing pixels with no cells.
#' 
#' @return The output is returned as a new \code{SpatialExperiment} object with 
#' \code{assay} slot containing the feature (cell types) x observations (pixels) matrix 
#' (dgCmatrix), \code{spatialCoords} slot containing spatial x,y coordinates of 
#' pixel centroids, and \code{colData} slot containing cell IDs of cells that 
#' were aggregated in each pixel.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom Matrix sparse.model.matrix
#' 
#' @export
#' 
rasterizeCellType <- function(input, col_name, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL, na.rm = FALSE) {
  ## extract cell type labels from SpatialExperiment
  cellTypes <- as.factor(colData(input)[,col_name])
  
  ## one-hot encode cell type labels as sparse matrix
  mat_ct <- t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
  rownames(mat_ct) <- levels(cellTypes)
  colnames(mat_ct) <- rownames(spatialCoords(input))
  
  ## rasterize
  out <- rasterizeSparseMatrix(mat_ct, spatialCoords(input), resolution = resolution, fun = fun, n_threads = 1, BPPARAM = BPPARAM)
  data_rast <- out$data_rast
  pos_rast <- out$pos_rast
  meta_rast <- out$meta_rast
  
  ## remove NA based on the na.rm argument
  if (na.rm) {
    na_cols <- colSums(is.na(data_rast)) != 0
    data_rast <- data_rast[,!na_cols]
    pos_rast <- pos_rast[!na_cols,]
    meta_rast <- meta_rast[!na_cols,]
  }
  
  ## construct a new SpatialExperiment object as output
  output <- SpatialExperiment::SpatialExperiment(
    assays = list(pixelval = data_rast),
    spatialCoords = pos_rast,
    colData = meta_rast
  )
  
  return(output)
}