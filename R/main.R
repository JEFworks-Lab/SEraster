#' rasterizeMatrix
#' 
#' @description Function to rasterize a given input matrix (both dense or sparse) based on a given position matrix.
#' 
#' @description This function assumes that inputs are provided as a \code{dgCmatrix} 
#' or \code{matrix} for data and \code{matrix} for position.
#' 
#' @param data \code{dgCmatrix} or \code{matrix}: Feature x observation matrix represented as a 
#' \code{dgCmatrix} or \code{matrix} object. Features can be genes or cell types. In the case of 
#' features being cell types, this matrix is assumed to be a sparse model matrix 
#' with rows as cell types and columns as cell IDs.
#' 
#' @param pos \code{matrix}: Spatial x,y coordinates of observations stored as a 
#' matrix array. Further, x,y coordinates are assumed to be stored in column 1 
#' and 2 of \code{spatialCoords}.
#' 
#' @param bbox \code{bbox}: Bounding box for rasterization defined by a numeric 
#' vector of length four, with xmin, ymin, xmax and ymax values.
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
#' x observation matrix as \code{dgCmatrix} if data was given as \code{dgCmatrix} 
#' and as \code{matrix} if data was given as \code{matrix}, spatial x,y coordinates of pixel 
#' centroids as \code{matrix}, and \code{data.frame} containing cell IDs of cells 
#' that were aggregated in each pixel.
#' 
#' @importFrom sf st_make_grid st_coordinates st_centroid st_as_sf st_intersects
#' @importFrom BiocParallel MulticoreParam bpstart bplapply bpstop
#' @importFrom Matrix rowMeans rowSums
#' 
#' @export
#' 
rasterizeMatrix <- function(data, pos, bbox, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL) {
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_threads)
  }
  BiocParallel::bpstart(BPPARAM)
  
  ## create grid for rasterization
  grid <- sf::st_make_grid(bbox, cellsize = resolution)
  
  ## extract position
  pos_pixel <- sf::st_coordinates(sf::st_centroid(grid))
  colnames(pos_pixel) <- c("x", "y")
  rownames(pos_pixel) <- paste0("pixel",seq_along(pos_pixel[,1]))
  
  ## convert pos to sf
  cells <- sf::st_as_sf(data.frame(pos), coords = c(1,2))
  ## get pixel ID for each cell
  ## since some cells are assigned to > 1 pixels if they are on the border, choose the 1st pixel id
  pixel_ids <- unlist(lapply(sf::st_intersects(cells, grid), function(sublist) sublist[1]))
  names(pixel_ids) <- rownames(pos)
  
  ## store aggregated subsetted matrix data for each pixel, store cell IDs and number of cells for each pixel into a data frame
  out <- unlist(BiocParallel::bplapply(sort(unique(pixel_ids)), function(id){
    ## get cell IDs for a particular pixel
    cell_ids <- names(pixel_ids[pixel_ids == id])
    
    ## subset feature observation matrix
    sub <- data[,cell_ids, drop = FALSE]
    ## aggregate cell counts to create pixel value
    if (fun == "mean") {
      pixel_val <- rowMeans(sub)
    } else if (fun == "sum") {
      pixel_val <- rowSums(sub)
    }
    
    ## store number of cells
    meta_rast <- data.frame(num_cell = length(cell_ids))
    ## store a list of cell IDs
    meta_rast$cellID_list <- list(cell_ids)
    
    if (is.matrix(data)) {
      return(list(pixel_val, meta_rast))
    } else {
      return(list(as(pixel_val, "CsparseMatrix"), meta_rast))
    }
  }), recursive = FALSE)
  
  ## stop parallel execution back-end
  BiocParallel::bpstop(BPPARAM)
  
  ## extract rasterized sparse matrix
  data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
  
  ## extract rasterized data frame
  meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
  
  ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
  rownames(data_rast) <- rownames(data)
  colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  
  ## subset rasterized pos
  pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast),]
  
  ## output
  output <- list("data_rast" = data_rast, "pos_rast" = pos_pixel, "meta_rast" = meta_rast)
}

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
#' @importFrom BiocParallel MulticoreParam bpstart bplapply bpstop
#' @importFrom Matrix rowMeans rowSums
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
  out <- unlist(BiocParallel::bplapply(sort(unique(pixel_ids)), function(id){
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
  
  ## stop parallel execution back-end
  BiocParallel::bpstop(BPPARAM)
  
  ## extract rasterized sparse matrix
  data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
  
  ## extract rasterized data frame
  meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
  
  ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
  rownames(data_rast) <- rownames(data)
  colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  
  ## subset rasterized pos
  pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast),]
  
  ## output
  output <- list("data_rast" = data_rast, "pos_rast" = pos_pixel, "meta_rast" = meta_rast)
}

#' rasterizeSparseMatrix2
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
#' @importFrom BiocParallel MulticoreParam bpstart bplapply bpstop bpoptions
#' @importFrom Matrix rowMeans rowSums
#' 
rasterizeSparseMatrix2 <- function(data, pos, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL) {
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_threads)
  }
  
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
  out <- unlist(BiocParallel::bplapply(sort(unique(pixel_ids)), function(id){
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
  }, BPPARAM = BPPARAM, BPOPTIONS = BiocParallel::bpoptions(packages = "Matrix")), recursive = FALSE)
  
  ## extract rasterized sparse matrix
  data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
  
  ## extract rasterized data frame
  meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
  
  ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
  rownames(data_rast) <- rownames(data)
  colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  
  ## subset rasterized pos
  pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast),]
  
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
#' @param input \code{SpatialExperiment} or \code{list}: Input data represented as a 
#' \code{SpatialExperiment} or \code{list} of \code{SpatialExperiment}. 
#' Each \code{SpatialExperiment} is assumed to have an \code{assay} slot containing feature (genes) x observation (cells)
#' matrix as \code{dgCmatrix} or \code{matrix} and a \code{spatialCoords} slot containing spatial x,y 
#' coordinates of observations as matrix array. Further, x,y coordinates are assumed 
#' to be stored in column 1 and 2 of \code{spatialCoords}.
#' 
#' @param assay_name \code{character}: Name of the assay slot of the input that 
#' you want to apply rasterization. If no argument is given, the first assay of the input 
#' would be rasterized. This argument is useful when you have both raw and normalized 
#' assays stored in the input, and you want to apply rasterization to the normalized assay. 
#' If the input is a \code{list}, assay_name is assumed to be present in all elements 
#' (\code{SpatialExperiment}) of the input.
#' 
#' @param resolution \code{integer}: Resolution or side length of each pixel. 
#' The unit of this parameter is assumed to be the same as the unit of spatial 
#' coordinates of the input data.
#' 
#' @param fun \code{character}: If "mean", pixel value for each pixel 
#' would be mean of gene expression for all cells within the pixel. If 
#' "sum", pixel value for each pixel would be sum of gene expression for all 
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
#' @return If the input was given as \code{SpatialExperiment}, the output is returned 
#' as a new \code{SpatialExperiment} object with \code{assay} slot containing the 
#' feature (genes) x observations (pixels) matrix (\code{dgCMatrix} or \code{matrix} 
#' depending on the input, see documentation for \code{rasterizeMatrix}), \code{spatialCoords} 
#' slot containing spatial x,y coordinates of pixel centroids, and \code{colData} 
#' slot containing cell IDs of cells that were aggregated in each pixel. If the input 
#' was provided as \code{list} of \code{SpatialExperiment}, the output is returned 
#' as a new \code{list} of \code{SpatialExperiment} containing information described 
#' above for corresponding \code{SpatialExperiment}. Further, \code{names(input)} 
#' is inherited in the output.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix colSums
#' 
#' @export
#' 
rasterizeGeneExpression <- function(input, assay_name = NULL, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- spatialCoords(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[i]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, xmin = min(pos[,1]), xmax = max(pos[,1]), ymin = min(pos[,2]), ymax = max(pos[,2])))
    }))
    bbox_common <- sf::st_bbox(c(
      xmin = floor(min(bbox_mat$xmin)-resolution/2), 
      xmax = ceiling(max(bbox_mat$xmax)+resolution/2), 
      ymin = floor(min(bbox_mat$ymin)-resolution/2), 
      ymax = ceiling(max(bbox_mat$ymax)+resolution/2)
    ))
    
    ## rasterize iteratively
    output_list <- lapply(seq_along(input), function(i) {
      ## get SpatialExperiment object of the given index
      spe <- input[[i]]
      
      if (is.null(assay_name)) {
        out <- SEraster::rasterizeMatrix(assay(spe), spatialCoords(spe), bbox = bbox_common, resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
      } else {
        stopifnot(is.character(assay_name))
        stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% assayNames(spe))
        out <- SEraster::rasterizeMatrix(assay(spe, assay_name), spatialCoords(spe), bbox = bbox_common, resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
      }
      data_rast <- out$data_rast
      pos_rast <- out$pos_rast
      meta_rast <- out$meta_rast
      
      ## construct a new SpatialExperiment object as output
      output <- SpatialExperiment::SpatialExperiment(
        assays = list(pixelval = data_rast),
        spatialCoords = pos_rast,
        colData = meta_rast
      )
      
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    
    ## return a list of SpatialExperiment object
    return(output_list)
    
  } else {
    ## create bbox
    pos <- spatialCoords(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## rasterize
    if (is.null(assay_name)) {
      out <- rasterizeMatrix(assay(input), spatialCoords(input), bbox = bbox, resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
    } else {
      stopifnot(is.character(assay_name))
      out <- rasterizeMatrix(assay(input, assay_name), spatialCoords(input), bbox = bbox, resolution = resolution, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM)
    }
    data_rast <- out$data_rast
    pos_rast <- out$pos_rast
    meta_rast <- out$meta_rast
    
    ## construct a new SpatialExperiment object as output
    output <- SpatialExperiment::SpatialExperiment(
      assays = list(pixelval = data_rast),
      spatialCoords = pos_rast,
      colData = meta_rast
    )
    
    ## return a SpatialExperiment object
    return(output)
  }
}

#' rasterizeCellType
#' 
#' @description Function to rasterize cell type labels in spatially-resolved 
#' transcriptomics data represented as SpatialExperiment class.
#'  
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} object.
#' 
#' @param input \code{SpatialExperiment} or \code{list}: Input data represented as a 
#' \code{SpatialExperiment} or \code{list} of \code{SpatialExperiment}. 
#' Each \code{SpatialExperiment} is assumed to have a \code{colData} slot containing cell type 
#' labels for observations as a data frame column and a \code{spatialCoords} slot 
#' containing spatial x,y coordinates of observations as matrix array. Further, 
#' x,y coordinates are assumed to be stored in column 1 and 2 of \code{spatialCoords}.
#' 
#' @param col_name \code{character}: Column name of the \code{colData} object 
#' containing cell type labels for observations. If the input is a \code{list}, 
#' col_name is assumed to be present in all elements (\code{SpatialExperiment}) of the input.
#' 
#' @param resolution \code{integer}: Resolution or side length of each pixel. 
#' The unit of this parameter is assumed to be the same as the unit of spatial 
#' coordinates of the input data.
#' 
#' @param fun \code{character}: If "mean", pixel value for each pixel 
#' would be the proportion of each cell type based on the one-hot-encoded cell type 
#' labels for all cells within the pixel. If "sum", pixel value for each pixel would 
#' be the number of cells of each cell type based on the one-hot-encoded cell type 
#' labels for all cells within the pixel.
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
#' @return If the input was given as \code{SpatialExperiment}, the output is returned 
#' as a new \code{SpatialExperiment} object with \code{assay} slot containing the 
#' feature (cell types) x observations (pixels) matrix (dgCmatrix), \code{spatialCoords} 
#' slot containing spatial x,y coordinates of pixel centroids, and \code{colData} 
#' slot containing cell IDs of cells that were aggregated in each pixel. If the input 
#' was provided as \code{list} of \code{SpatialExperiment}, the output is returned 
#' as a new \code{list} of \code{SpatialExperiment} containing information described 
#' above for corresponding \code{SpatialExperiment}. Further, \code{names(input)} 
#' is inherited in the output.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom Matrix sparse.model.matrix t colSums
#' 
#' @export
#' 
rasterizeCellType <- function(input, col_name, resolution = 100, fun = "mean", n_threads = 1, BPPARAM = NULL) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- spatialCoords(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[i]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, xmin = min(pos[,1]), xmax = max(pos[,1]), ymin = min(pos[,2]), ymax = max(pos[,2])))
    }))
    bbox_common <- sf::st_bbox(c(
      xmin = floor(min(bbox_mat$xmin)-resolution/2), 
      xmax = ceiling(max(bbox_mat$xmax)+resolution/2), 
      ymin = floor(min(bbox_mat$ymin)-resolution/2), 
      ymax = ceiling(max(bbox_mat$ymax)+resolution/2)
    ))
    
    ## rasterize iteratively
    output_list <- lapply(seq_along(input), function(i) {
      ## get SpatialExperiment object of the given index
      spe <- input[[i]]
      
      ## extract cell type labels from SpatialExperiment
      stopifnot(is.character(col_name))
      stopifnot("col_name does not exist in the input SpatialExperiment object"= col_name %in% colnames(colData(spe)))
      cellTypes <- as.factor(colData(spe)[,col_name])
      
      ## one-hot encode cell type labels as sparse matrix
      mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
      rownames(mat_ct) <- levels(cellTypes)
      colnames(mat_ct) <- rownames(spatialCoords(spe))
      
      ## rasterize
      out <- rasterizeMatrix(mat_ct, spatialCoords(spe), bbox_common, resolution = resolution, fun = fun, n_threads = 1, BPPARAM = BPPARAM)
      data_rast <- out$data_rast
      pos_rast <- out$pos_rast
      meta_rast <- out$meta_rast
      
      ## construct a new SpatialExperiment object as output
      output <- SpatialExperiment::SpatialExperiment(
        assays = list(pixelval = data_rast),
        spatialCoords = pos_rast,
        colData = meta_rast
      )
      
      return(output)
    })
    
    if (!is.null(names(input))) {
      names(output_list) <- names(input)
    }
    
    ## return a list of SpatialExperiment object
    return(output_list)
    
  } else {
    ## create bbox
    pos <- spatialCoords(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## extract cell type labels from SpatialExperiment
    stopifnot(is.character(col_name))
    stopifnot("col_name does not exist in the input SpatialExperiment object"= col_name %in% colnames(colData(spe)))
    cellTypes <- as.factor(colData(input)[,col_name])
    
    ## one-hot encode cell type labels as sparse matrix
    mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
    rownames(mat_ct) <- levels(cellTypes)
    colnames(mat_ct) <- rownames(spatialCoords(input))
    
    ## rasterize
    out <- rasterizeMatrix(mat_ct, spatialCoords(input), bbox, resolution = resolution, fun = fun, n_threads = 1, BPPARAM = BPPARAM)
    data_rast <- out$data_rast
    pos_rast <- out$pos_rast
    meta_rast <- out$meta_rast
    
    ## construct a new SpatialExperiment object as output
    output <- SpatialExperiment::SpatialExperiment(
      assays = list(pixelval = data_rast),
      spatialCoords = pos_rast,
      colData = meta_rast
    )
    
    ## return a SpatialExperiment object
    return(output)
  }
}

#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix colSums colMeans
#' @importFrom ggplot2 ggplot theme ggtitle
#' 
#' @export
plotRaster <- function(input, assay_name = NULL, feature_name = "sum", factor_levels = NULL, showLegend = TRUE, plotTitle = NULL, showMinimal = TRUE, ...) {
  ## get the indicated assay slot (features-by-observations matrix)
  if (is.null(assay_name)) {
    mat <- assay(input)
  } else {
    stopifnot(is.character(assay_name))
    stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% assayNames(input))
    mat <- assay(input, assay_name)
  }
  
  ## create data.frame for plotting
  if (feature_name == "sum") {
    df <- data.frame(x = SpatialExperiment::spatialCoords(input)[,1], y = SpatialExperiment::spatialCoords(input)[,2], fill = colSums(mat))
  } else if (feature_name == "mean") {
    df <- data.frame(x = SpatialExperiment::spatialCoords(input)[,1], y = SpatialExperiment::spatialCoords(input)[,2], fill = colMeans(mat))
  } else {
    stopifnot(is.character(feature_name))
    stopifnot("feature_name does not exist in the input SpatialExperiment object's assay slot" = feature_name %in% rownames(mat))
    df <- data.frame(x = SpatialExperiment::spatialCoords(input)[,1], y = SpatialExperiment::spatialCoords(input)[,2], fill = mat[feature_name,])
  }
  
  ## compute resolution
  resolution <- diff(spatialCoords(input))[1,1]
  
  ## change object class of fill if plotting categorical variables
  if (is.null(factor_levels)) {
    plt <- ggplot2::ggplot(df, aes(x = x, y = y, fill = fill)) +
      coord_fixed() +
      geom_tile(width = resolution, height = resolution) +
      scale_fill_viridis_c(...) +
      theme_bw() +
      theme(panel.grid = element_blank())
  } else {
    df$fill <- factor(df$fill, levels = factor_levels)
    plt <- ggplot2::ggplot(df, aes(x = x, y = y, fill = fill)) +
      coord_fixed() +
      geom_tile(width = resolution, height = resolution) +
      scale_fill_viridis_d(...) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
  
  if (!showLegend) {
    plt <- plt + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(plotTitle)) {
    plt <- plt + ggplot2::ggtitle(plotTitle)
  }
  
  if (showMinimal) {
    plt <- plt + ggplot2::theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  }
  
  return(plt)
}