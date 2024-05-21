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
#' @param bbox \code{bbox} or \code{numeric}: Bounding box for rasterization defined 
#' by a bbox class object (as created by sf::st_bbox) or a numeric vector of length four, with xmin, ymin, xmax 
#' and ymax values. Values in a numeric vector are assumed to be in the order of xmin, 
#' ymin, xmax, and ymax.
#' 
#' @param resolution \code{integer} or \code{double}: Resolution refers to the side 
#' length of each pixel for square pixels and the distance between opposite edges 
#' of each pixel for hexagonal pixels. The unit of this parameter is assumed to 
#' be the same as the unit of spatial coordinates of the input data.
#' 
#' @param square \code{logical}: If TRUE (default), rasterize into square pixels. If FALSE, 
#' rasterize into hexagonal pixels.
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
#' @param verbose \code{logical}: Whether to display verbose output or warning. Default is FALSE 
#' 
#' @return The output is returned as a \code{list} containing rasterized feature 
#' x observation matrix as \code{dgCmatrix} if data was given as \code{dgCmatrix} 
#' and as \code{matrix} if data was given as \code{matrix}, spatial x,y coordinates of pixel 
#' centroids as \code{matrix}, and \code{data.frame} containing meta data for pixels 
#' (number of cells that were aggregated in each pixel, cell IDs of cells that were 
#' aggregated in each pixel, pixel type based on the \code{square} argument, pixel 
#' resolution based on the \code{resolution} argument, pixel geometry as \code{sfc_POLYGON}).
#' 
#' @importFrom sf st_make_grid st_coordinates st_centroid st_as_sf st_intersects
#' @importFrom BiocParallel MulticoreParam bpstart bplapply bpstop
#' @importFrom Matrix rowMeans rowSums
#' @importFrom methods as
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(sf)
#' 
#' data("merfish_mousePOA")
#' # extract features-by-cells matrix, spatial coordinates from the SpatialExperiment object
#' data <- assay(merfish_mousePOA)
#' pos <- spatialCoords(merfish_mousePOA)
#' # compute bounding box
#' resolution <- 100
#' bbox <- st_bbox(c(
#'   xmin = floor(min(pos[,1])-resolution/2), 
#'   xmax = ceiling(max(pos[,1])+resolution/2), 
#'   ymin = floor(min(pos[,2])-resolution/2), 
#'   ymax = ceiling(max(pos[,2])+resolution/2)
#' ))
#' # rasterize with mean as the aggregation function
#' out_mean <- rasterizeMatrix(data, pos, bbox, resolution = resolution, fun = "mean")
#' # rasterize with sum as the aggregation function
#' out_sum <- rasterizeMatrix(data, pos, bbox, resolution = resolution, fun = "sum")
#' # rasterize with user-defined resolution and hexagonal pixels
#' # in this case, you need to update the bbox as well
#' resolution <- 200
#' bbox <- st_bbox(c(
#'   xmin = floor(min(pos[,1])-resolution/2), 
#'   xmax = ceiling(max(pos[,1])+resolution/2), 
#'   ymin = floor(min(pos[,2])-resolution/2), 
#'   ymax = ceiling(max(pos[,2])+resolution/2)
#' ))
#' out_hex <- rasterizeMatrix(data, pos, bbox, resolution = resolution, square = FALSE, fun = "mean")
#' 
rasterizeMatrix <- function(data, pos, bbox, resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = TRUE) {
  ## set up parallel execution back-end with BiocParallel
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = n_threads)
  }
  # BiocParallel::bpstart(BPPARAM)
  
  ## check to see if data and pos are have the same number of cells
  stopifnot("data and pos must have same number of cells" = ncol(data)==nrow(pos))
  
  ## check if cellnames in data and pos are identical if they are provided in colnames and rownames respectively
  stopifnot((is.null(colnames(data)) | is.null(rownames(pos))) | identical(colnames(data), rownames(pos)))|
  
  ## check to see if bbox input is of class numeric or bbox, if numeric, convert to st_bbox
  stopifnot("bbox must be of class 'numeric' or 'bbox' (as generated by sf::st_bbox)" = class(bbox) %in% c('numeric', 'bbox'))
  if (class(bbox)=='numeric') {
    stopifnot("Bounding box for rasterization must be defined by a numeric vector of length four, with xmin, ymin, xmax and ymax values or object of class bbox (as generated by sf::st_bbox)" = length(bbox)==4)
    bbox <- sf::st_bbox(c(xmin = bbox[1], ymin = bbox[2], xmax = bbox[3], ymax = bbox[4]))
  } 
  ## create grid for rasterization
  grid <- sf::st_make_grid(bbox, cellsize = resolution, square = square)
  
  if (verbose) {
    if (resolution >= diff(range(pos[,1])) | resolution >= diff(range(pos[,2]))) {
      warning("Resolution is equal to or larger than the input data's x or y dimensions. Consider using lower resolutions.")
    }
  }
  
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
  }, BPPARAM = BPPARAM), recursive = FALSE)
  
  ## stop parallel execution back-end
  # BiocParallel::bpstop(BPPARAM)
  
  ## extract rasterized sparse matrix
  data_rast <- do.call(cbind, out[seq(1,length(out),by=2)])
  
  ## extract rasterized data frame
  meta_rast <- do.call(rbind, out[seq(1,length(out),by=2)+1])
  
  ## set rownames/colnames for rasterized sparse matrix and rasterized data frame
  rownames(data_rast) <- rownames(data)
  colnames(data_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  rownames(meta_rast) <- paste0("pixel", sort(unique(pixel_ids)))
  
  ## subset rasterized pos
  pos_pixel <- pos_pixel[rownames(pos_pixel) %in% colnames(data_rast), , drop = FALSE]
  
  ## add pixel information to data frame
  # pixel type
  if (square) {
    meta_rast$type <- "square"
  } else {
    meta_rast$type <- "hexagon"
  }
  # pixel resolution
  meta_rast$resolution <- resolution
  # pixel geometry
  df_sf <- sf::st_sf(geometry = grid, row.names = paste0("pixel",seq_along(grid)))
  df_sf <- df_sf[rownames(df_sf) %in% rownames(pos_pixel),]
  meta_rast$geometry <- df_sf$geometry
  
  ## output
  output <- list("data_rast" = data_rast, "pos_rast" = pos_pixel, "meta_rast" = meta_rast)
}


#' rasterizeGeneExpression
#' 
#' @description Function to rasterize feature x observation matrix in spatially-resolved 
#' omics data represented as SpatialExperiment class.
#'  
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.
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
#' @param resolution \code{integer} or \code{double}: Resolution refers to the side 
#' length of each pixel for square pixels and the distance between opposite edges 
#' of each pixel for hexagonal pixels. The unit of this parameter is assumed to 
#' be the same as the unit of spatial coordinates of the input data.
#' 
#' @param square \code{logical}: If TRUE (default), rasterize into square pixels. If FALSE, 
#' rasterize into hexagonal pixels.
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
#' @param verbose \code{logical}: Whether to display verbose output or warning. Default is FALSE 
#' 
#' @return If the input was given as \code{SpatialExperiment}, the output is returned 
#' as a new \code{SpatialExperiment} object with \code{assay} slot containing the 
#' feature (genes) x observations (pixels) matrix (\code{dgCMatrix} or \code{matrix} 
#' depending on the input, see documentation for \link{rasterizeMatrix}), \code{spatialCoords} 
#' slot containing spatial x,y coordinates of pixel centroids, and \code{colData} slot 
#' containing meta data for pixels (number of cells that were aggregated in each pixel, 
#' cell IDs of cells that were aggregated in each pixel, pixel type based on the 
#' \code{square} argument, pixel resolution based on the \code{resolution} argument, 
#' pixel geometry as \code{sfc_POLYGON}). If the input was provided as \code{list} 
#' of \code{SpatialExperiment}, the output is returned as a new \code{list} of 
#' \code{SpatialExperiment} containing information described above for corresponding 
#' \code{SpatialExperiment}. Further, \code{names(input)} is inherited in the output.

#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom Matrix colSums
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' 
#' data("merfish_mousePOA")
#' 
#' # check assay names for this particular SpatialExperiment object (should be "volnorm")
#' assayNames(merfish_mousePOA)
#' 
#' # rasterize a single SpatialExperiment object
#' # make sure to specify the assay_name argument when the input SpatialExperiment 
#' # object has multiple assay names (assay_name is used here as an example)
#' out <- rasterizeGeneExpression(merfish_mousePOA, assay_name = "volnorm", fun = "mean")
#' 
#' # rasterize a single SpatialExperiment object with user-defined resolution and hexagonal pixels
#' out <- rasterizeGeneExpression(merfish_mousePOA, assay_name = "volnorm", resolution = 200, square = FALSE, fun = "mean")
#' 
#' # rasterize a list of SpatialExperiment objects (in this case, permutated datasets 
#' # with 3 different rotations)
#' spe_list <- permutateByRotation(merfish_mousePOA, n_perm = 3)
#' out_list <- rasterizeGeneExpression(spe_list, assay_name = "volnorm", resolution = 100, square = TRUE, fun = "mean")
#' 
rasterizeGeneExpression <- function(input, assay_name = NULL, resolution = 100, square = TRUE, fun = "mean", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- SpatialExperiment::spatialCoords(input[[i]])
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
        out <- SEraster::rasterizeMatrix(SummarizedExperiment::assay(spe), SpatialExperiment::spatialCoords(spe), bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
      } else {
        stopifnot(is.character(assay_name))
        stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% SummarizedExperiment::assayNames(spe))
        out <- SEraster::rasterizeMatrix(SummarizedExperiment::assay(spe, assay_name), SpatialExperiment::spatialCoords(spe), bbox = bbox_common, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
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
    pos <- SpatialExperiment::spatialCoords(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## rasterize
    if (is.null(assay_name)) {
      out <- rasterizeMatrix(SummarizedExperiment::assay(input), SpatialExperiment::spatialCoords(input), bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
    } else {
      stopifnot(is.character(assay_name))
      stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% SummarizedExperiment::assayNames(input))
      out <- rasterizeMatrix(SummarizedExperiment::assay(input, assay_name), SpatialExperiment::spatialCoords(input), bbox = bbox, resolution = resolution, square = square, fun = fun, n_threads = n_threads, BPPARAM = BPPARAM, verbose = verbose)
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
#' omics data represented as SpatialExperiment class.
#'
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.
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
#' @param resolution \code{integer} or \code{double}: Resolution refers to the side 
#' length of each pixel for square pixels and the distance between opposite edges 
#' of each pixel for hexagonal pixels. The unit of this parameter is assumed to 
#' be the same as the unit of spatial coordinates of the input data.
#' 
#' @param square \code{logical}: If TRUE (default), rasterize into square pixels. If FALSE, 
#' rasterize into hexagonal pixels.
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
#' @param verbose \code{logical}: Whether to display verbose output or warning. Default is FALSE 
#' 
#' @return If the input was given as \code{SpatialExperiment}, the output is returned 
#' as a new \code{SpatialExperiment} object with \code{assay} slot containing the 
#' feature (cell types) x observations (pixels) matrix (dgCmatrix), \code{spatialCoords} 
#' slot containing spatial x,y coordinates of pixel centroids, and \code{colData} slot 
#' containing meta data for pixels (number of cells that were aggregated in each pixel, 
#' cell IDs of cells that were aggregated in each pixel, pixel type based on the 
#' \code{square} argument, pixel resolution based on the \code{resolution} argument, 
#' pixel geometry as \code{sfc_POLYGON}). If the input was provided as \code{list} 
#' of \code{SpatialExperiment}, the output is returned as a new \code{list} of 
#' \code{SpatialExperiment} containing information described above for corresponding 
#' \code{SpatialExperiment}. Further, \code{names(input)} is inherited in the output.
#' 
#' @importFrom SpatialExperiment spatialCoords SpatialExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom Matrix sparse.model.matrix t colSums
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' 
#' data("merfish_mousePOA")
#' 
#' # check assay names for this particular SpatialExperiment object (you can see 
#' # that cell-type labels are stored in the "celltype" column)
#' head(colData(merfish_mousePOA))
#' 
#' # rasterize a single SpatialExperiment object
#' # make sure to specify the col_name argument
#' out <- rasterizeCellType(merfish_mousePOA, col_name = "celltype", fun = "sum")
#' 
#' # rasterize a single SpatialExperiment object with user-defined resolution and hexagonal pixels
#' out <- rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 200, square = FALSE, fun = "sum")
#' 
#' # rasterize a list of SpatialExperiment objects (in this case, permutated datasets 
#' # with 3 different rotations)
#' spe_list <- permutateByRotation(merfish_mousePOA, n_perm = 3)
#' out_list <- rasterizeCellType(spe_list, col_name = "celltype", resolution = 100, square = TRUE, fun = "sum")
#' 
rasterizeCellType <- function(input, col_name, resolution = 100, square = TRUE, fun = "sum", n_threads = 1, BPPARAM = NULL, verbose = FALSE) {
  if (is.list(input)) {
    ## create a common bbox
    bbox_mat <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- SpatialExperiment::spatialCoords(input[[i]])
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
      colnames(mat_ct) <- rownames(SpatialExperiment::spatialCoords(spe))
      
      ## rasterize
      out <- rasterizeMatrix(mat_ct, SpatialExperiment::spatialCoords(spe), bbox_common, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
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
    pos <- SpatialExperiment::spatialCoords(input)
    bbox <- sf::st_bbox(c(
      xmin = floor(min(pos[,1])-resolution/2), 
      xmax = ceiling(max(pos[,1])+resolution/2), 
      ymin = floor(min(pos[,2])-resolution/2), 
      ymax = ceiling(max(pos[,2])+resolution/2)
    ))
    
    ## extract cell type labels from SpatialExperiment
    stopifnot(is.character(col_name))
    stopifnot("col_name does not exist in the input SpatialExperiment object"= col_name %in% colnames(colData(input)))
    cellTypes <- as.factor(colData(input)[,col_name])
    
    ## one-hot encode cell type labels as sparse matrix
    mat_ct <- Matrix::t(Matrix::sparse.model.matrix(~ 0 + cellTypes))
    rownames(mat_ct) <- levels(cellTypes)
    colnames(mat_ct) <- rownames(SpatialExperiment::spatialCoords(input))
    
    ## rasterize
    out <- rasterizeMatrix(mat_ct, SpatialExperiment::spatialCoords(input), bbox, resolution = resolution, square = square, fun = fun, n_threads = 1, BPPARAM = BPPARAM, verbose = verbose)
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


#' permutateByRotation
#' 
#' @description Function to permutate a given input SpatialExperiment object(s) 
#' by rotating the x,y coordinates around the midrange point.
#' 
#' @description This function assumes that the input is provided as a \code{SpatialExperiment} 
#' object or a \code{list} of \code{SpatialExperiment} objects.
#'
#' @description When the input is a \code{list} of \code{SpatialExperiment} objects, 
#' all \code{SpatialExperiment} objects will be rotated around a common midrange 
#' point computed based on combined x,y coordinates.
#' 
#' @param input \code{SpatialExperiment} or \code{list}: Input data represented as a 
#' \code{SpatialExperiment} or \code{list} of \code{SpatialExperiment}. 
#' Each \code{SpatialExperiment} is assumed to have an \code{assay} slot containing feature (genes) x observation (cells)
#' matrix as \code{dgCmatrix} or \code{matrix} and a \code{spatialCoords} slot containing spatial x,y 
#' coordinates of observations as matrix array. Further, x,y coordinates are assumed 
#' to be stored in column 1 and 2 of \code{spatialCoords}, and column names of \code{spatialCoords} 
#' are assumed to be "x" and "y", respectively.
#' 
#' @param n_perm \code{integer}: Number of permutations. Default = 1. This number is used to compute the angles at which the input is rotated at.
#' 
#' @param verbose \code{logical}: Whether to display verbose output or warning. Default is FALSE. 
#' 
#' @return If the input was given as \code{SpatialExperiment}, the output is returned 
#' as a \code{list} of \code{n_perm} \code{SpatialExperiment} objects. Each \code{SpatialExperiment} 
#' object has an updated \code{spatialCoords} slot containing the spatial x,y coordinates 
#' rotated at a corresponding angle. \code{assay} and \code{colData} slots are inherited.
#' Further, \code{names()} of the output indicates the angles at which the input 
#' is rotated at. If the input was given as \code{list} of \code{SpatialExperiment}, 
#' the output is returned as a new \code{list} of \code{length(input)} * \code{n_perm} 
#' \code{SpatialExperiment} objects. Each \code{SpatialExperiment} object has an 
#' updated \code{spatialCoords} slot containing the spatial x,y coordinates rotated 
#' at a corresponding angle. \code{assay} and \code{colData} slots are inherited.
#' Further, \code{names()} of the output indicates the dataset names from \code{names(input)} 
#' and the angles at which the input is rotated at.
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays colData
#' @importFrom rearrr rotate_2d midrange
#' 
#' @export
#' 
#' @examples
#' data("merfish_mousePOA")
#' 
#' # create a list of 3 permutated datasets rotated at 0 (original), 120, and 240 degrees
#' # this output can directly be fed into rasterizeGeneExpression or rasterizeCellType 
#' # functions to rasterize all 3 permutations at once with the same pixel coordinates
#' spe_list <- permutateByRotation(merfish_mousePOA, n_perm = 3)
#' 
#' # create a list of 5 permutated datasets rotated at 0 (original), 72, 144, 216, 288 degrees
#' spe_list <- permutateByRotation(merfish_mousePOA, n_perm = 5)
#' 
permutateByRotation <- function(input, n_perm = 1, verbose = FALSE) {
  ## compute rotation degrees based on the required number of permutation
  angles <- seq(0, 360, by = 360/n_perm)[1:n_perm]
  
  if (verbose) {
    message(paste0("Number of permutations = ", n_perm))
    message(paste0("Angles used for rotations are ", paste(angles, collapse = ", "), " degrees"))
  }
  
  if (is.list(input)) {
    ## combine all x,y coordinates
    pos_comb <- do.call(rbind, lapply(seq_along(input), function(i) {
      pos <- SpatialExperiment::spatialCoords(input[[i]])
      if (!is.null(names(input))) {
        dataset <- names(input)[i]
      } else {
        dataset <- i
      }
      return(data.frame(dataset = dataset, x = pos[,1], y = pos[,2]))
    }))
    ## find the midrange point across combined x,y coordinates
    midrange_pt <- rearrr::midrange(pos_comb, cols = c("x", "y"))
    
    if (verbose) {
      message(paste0("Rotating all datasets around (x, y) = (", midrange_pt$x, ", ", midrange_pt$y, ")."))
    }
    
    output <- unlist(lapply(input, function(spe) {
      stopifnot("input must be a SpatialExperiment object or a list of SpatialExperiment objects"=class(spe) == "SpatialExperiment")
      
      ## get original x,y coordinates
      pos_orig <- data.frame(spatialCoords(spe))
      stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))
      
      output2 <- lapply(angles, function(angle) {
        ## rotate around the midrange point
        pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin = as.numeric(midrange_pt), overwrite = TRUE)
        
        pos_rotated <- as.matrix(pos_rotated[,c("x_rotated", "y_rotated")])
        colnames(pos_rotated) <- c("x", "y")
        rownames(pos_rotated) <- rownames(pos_orig)
        
        ## update SpatialExperiment object
        spe_rotated <- SpatialExperiment::SpatialExperiment(
          assays = SummarizedExperiment::assays(spe),
          spatialCoords = pos_rotated,
          colData = SummarizedExperiment::colData(spe)
        )
        return(spe_rotated)
      })
      return(output2)
    }))
    
    ## assign list names
    if (!is.null(names(input))) {
      names(output) <- paste0(rep(names(input), each = length(angles)), "_rotated_", angles)
    }
    
    ## return a list of SpatialExperiment objects
    return(output)
    
  } else {
    stopifnot("input must be a SpatialExperiment object or a list of SpatialExperiment objects"=class(input) == "SpatialExperiment")
    
    ## get original x,y coordinates
    pos_orig <- data.frame(spatialCoords(input))
    stopifnot("Column 1 and 2 of the spatialCoords slot should be named x and y, respectively. Please change column names accordingly."=colnames(pos_orig)[1:2] == c("x", "y"))
    
    output  <- lapply(angles, function(angle) {
      ## rotate around the midrange point
      pos_rotated <- rearrr::rotate_2d(data = pos_orig, degrees = angle, x_col = "x", y_col = "y", origin_fn = rearrr::midrange, overwrite = TRUE)
      
      pos_rotated <- as.matrix(pos_rotated[,c("x_rotated", "y_rotated")])
      colnames(pos_rotated) <- c("x", "y")
      rownames(pos_rotated) <- rownames(pos_orig)
      
      ## update SpatialExperiment object
      spe_rotated <- SpatialExperiment::SpatialExperiment(
        assays = SummarizedExperiment::assays(input),
        spatialCoords = pos_rotated,
        colData = SummarizedExperiment::colData(input)
      )
      return(spe_rotated)
    })
    
    ## assign list names
    names(output) <- paste0("rotated_", angles)
    
    ## return a list of SpatialExperiment objects
    return(output)
  }
}


#' plotRaster
#' 
#' @description Function based on \code{ggplot2::geom_tile} to visualize a rasterized spatial omics dataset represented as a \code{SpatialExperiment} object.
#' 
#' @param input \code{SpatialExperiment}: Input data represented as a 
#' \code{SpatialExperiment}. The given \code{SpatialExperiment} is assumed to have 
#' an \code{assay} slot containing a features-by-observations matrix as \code{dgCmatrix} 
#' or \code{matrix} and a \code{colData} slot containing \code{sfc_POLYGON} geometry 
#' of pixels. The features-by-observations matrix is assumed to have either genes 
#' or cell types as features and pixels as observations.
#' 
#' @param assay_name \code{character}: Name of the assay slot of the input that 
#' you want to visualize. If no argument is given, the first assay of the input 
#' would be visualized. This argument is useful when you have multiple assays stored 
#' in the input, and you want to visualize a specific assay. Default is NULL.
#' 
#' @param feature_name \code{character}: Name of the feature in the input that you 
#' want to visualize. This argument is useful when you want to specify a feature you 
#' want to visualize. You can also use "sum" to visualize sum of all feature 
#' values per observation or "mean" to visualize mean of all feature values per observation. 
#' Default is "sum".
#' 
#' @param factor_levels \code{character} or \code{numeric} or \code{factor}: An 
#' optional vector to convert and plot the input data as \code{factor}. This argument 
#' is useful if you want to plot categorical/ordinal variables, such as binarized 
#' occurrence of a specific cell type. \code{factor_levels} is fed into \code{levels} 
#' argument of the \code{factor} function in base R. Default is NULL.
#' 
#' @param showLegend \code{logical}: Boolean to show the legend. Default is TRUE.
#' 
#' @param plotTitle \code{character}: An optional argument to add a title to the 
#' resulting plot. Default is NULL.
#' 
#' @param showAxis \code{logical}: Boolean to show axis title, texts, and ticks. 
#' Default is FALSE.
#' 
#' @param ... Additional parameters to pass to \code{ggplot2::scale_fill_viridis_c} 
#' if no argument is provided to \code{factor_levels} or \code{ggplot2::scale_fill_viridis_d} 
#' if a vector is provided to \code{factor_levels}. If you wish to use other color 
#' maps, we recommend overriding the resulting \code{ggplot} object.
#' 
#' @return The output is returned as a \code{ggplot} object, where the input is 
#' visualized as \code{ggplot2::geom_sf}. Each pixel is plotted based on \code{sfc_POLYGON} 
#' geometry stored in the \code{colData} slot. Coloring of pixel represent the corresponding 
#' values of summarized (sum or mean) or specific feature (e.g. rasterized gene expression) 
#' per observation (pixel).
#' 
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix colSums colMeans
#' @importFrom sf st_sf
#' @importFrom ggplot2 ggplot aes coord_fixed geom_sf scale_fill_viridis_c scale_fill_viridis_d theme_bw theme ggtitle element_blank
#' 
#' @export
#' 
#' @examples
#' data("merfish_mousePOA")
#' 
#' # rasterize gene expression
#' out <- rasterizeGeneExpression(merfish_mousePOA, assay_name = "volnorm", fun = "mean")
#' 
#' # plot total rasterized gene expression per pixel (there is only one assay_name 
#' # in out and default for feature_name argument is "sum"; therefore, these arguments 
#' # are not specified)
#' plotRaster(out, name = "total rasterized gexp")
#' 
#' # plot rasterized expression of a specific gene/feature per pixel
#' plotRaster(out, feature_name = "Esr1", name = "Esr1")
#' 
#' # rasterize cell-type labels with user-defined resolution and hexagonal pixels
#' out <- rasterizeCellType(merfish_mousePOA, col_name = "celltype", resolution = 50, square = FALSE, fun = "sum")
#' 
#' # plot total cell counts per pixel (there is only one assay_name in out and default 
#' # for feature_name argument is "sum"; therefore, these arguments are not specified)
#' # here, let's use additional parameters for ggplot2::scale_fill_viridis_c so 
#' # that it would have a different color scheme from gene expression plots
#' plotRaster(out, name = "total cell counts", option = "inferno")
#' 
#' # plot specific cell type's cell counts per pixel
#' plotRaster(out, feature_name = "Inhibitory", name = "Inhibitory neuron counts", option = "inferno")
#' 
plotRaster <- function(input, assay_name = NULL, feature_name = "sum", factor_levels = NULL, showLegend = TRUE, plotTitle = NULL, showAxis = FALSE, ...) {
  ## get the indicated assay slot (features-by-observations matrix)
  if (is.null(assay_name)) {
    mat <- SummarizedExperiment::assay(input)
  } else {
    stopifnot(is.character(assay_name))
    stopifnot("assay_name does not exist in the input SpatialExperiment object"= assay_name %in% SummarizedExperiment::assayNames(input))
    mat <- SummarizedExperiment::assay(input, assay_name)
  }
  
  ## create data.frame for plotting
  # create sf data.frame with geometry (sfc_POLYGON) from colData
  df_sf <- sf::st_sf(geometry = colData(input)$geometry, row.names = rownames(colData(input)))
  # add pixel values
  if (feature_name == "sum") {
    df_sf <- cbind(df_sf, fill_var = colSums(mat))
  } else if (feature_name == "mean") {
    df_sf <- cbind(df_sf, fill_var = colMeans(mat))
  } else {
    stopifnot(is.character(feature_name))
    stopifnot("feature_name does not exist in the input SpatialExperiment object's assay slot" = feature_name %in% rownames(mat))
    df_sf <- cbind(df_sf, fill_var = mat[feature_name,])
  }
  
  if (is.null(factor_levels)) {
    plt <- ggplot2::ggplot() +
      ggplot2::coord_fixed() +
      ggplot2::geom_sf(data = df_sf, ggplot2::aes(fill = fill_var)) +
      ggplot2::scale_fill_viridis_c(...) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  } else {
    df_sf$fill_var <- factor(df_sf$fill_var, levels = factor_levels)
    plt <- ggplot2::ggplot() +
      ggplot2::coord_fixed() +
      ggplot2::geom_sf(data = df_sf, ggplot2::aes(fill = fill_var)) +
      ggplot2::scale_fill_viridis_d(...) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  
  if (!showLegend) {
    plt <- plt + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.null(plotTitle)) {
    plt <- plt + ggplot2::ggtitle(plotTitle)
  }
  
  if (!showAxis) {
    plt <- plt + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  }
  
  return(plt)
}
