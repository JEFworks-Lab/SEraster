library(SpatialExperiment)
library(Matrix)
library(sf)

data("merfish_mousePOA")

resolution <- 100

data <- assay(merfish_mousePOA)
pos <- spatialCoords(merfish_mousePOA)
bbox <- sf::st_bbox(c(
  xmin = floor(min(pos[,1])-resolution/2), 
  xmax = ceiling(max(pos[,1])+resolution/2), 
  ymin = floor(min(pos[,2])-resolution/2), 
  ymax = ceiling(max(pos[,2])+resolution/2)
))

out_rasterizeMatrix <- rasterizeMatrix(data, pos, bbox)

test_that("Correct structure of rasterizeMatrix output (list)", {
  # output is a list
  expect_type(out_rasterizeMatrix, "list")
  # output has correct element names
  expect_equal(names(out_rasterizeMatrix), c("data_rast", "pos_rast", "meta_rast"))
  # output elements have the same number of pixels
  expect_equal(ncol(out_rasterizeMatrix$data_rast), nrow(out_rasterizeMatrix$pos_rast))
  expect_equal(ncol(out_rasterizeMatrix$data_rast), nrow(out_rasterizeMatrix$meta_rast))
  expect_equal(nrow(out_rasterizeMatrix$pos_rast), nrow(out_rasterizeMatrix$meta_rast))
  # output element (data_rast) retained the same number of genes
  expect_equal(nrow(out_rasterizeMatrix$data_rast), nrow(merfish_mousePOA))
  # output element (pos_rast) has x,y coordinates
  expect_equal(colnames(out_rasterizeMatrix$pos_rast), c("x","y"))
  # output number of pixels ≤ number of cells in the input
  
})
