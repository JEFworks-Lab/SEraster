#' Preprocessed MERFISH dataset of the mouse preoptic area for a bregma -0.29 slice 
#' from a female naive animal (Animal ID = 1, Animal Sex = "Female", Behavior = 
#' "Naive", Bregma = "-0.29").
#'
#' @format \code{SpatialExperiment} object where \code{assay} slot contains genes-by-cells 
#' matrix with preprocessed gene expression (total RNA counts per cell divided by 
#' cell volume and scaled by 1000) as \code{dgCMatrix}, \code{spatialCoords} slot 
#' contains x,y coordinates of cells, and \code{colData} slot contains bregma, 
#' cell type, and neuron type meta data.
#'
#' @source \url{https://www.science.org/doi/10.1126/science.aau5324}
#' 
#' @usage data("merfish_mousePOA")
#' 
#' @return \code{SpatialExperiment} object for the preprocessed MERFISH dataset 
#' of the mouse preoptic area for a bregma -0.29 slice from a female naive animal 
#' (Animal ID = 1, Animal Sex = "Female", Behavior = "Naive", Bregma = "-0.29").
"merfish_mousePOA"