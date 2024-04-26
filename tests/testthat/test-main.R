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

out_rasterizeGeneExpression_spe <- rasterizeGeneExpression(merfish_mousePOA)

out_rasterizeCellType_spe <- rasterizeCellType(merfish_mousePOA, col_name = "celltype")

out_permutateByRotation <- permutateByRotation(merfish_mousePOA, n_perm = 10)

out_rasterizeGeneExpression_list <- rasterizeGeneExpression(out_permutateByRotation)

out_rasterizeCellType_list <- rasterizeCellType(out_permutateByRotation, col_name = "celltype")

out_plotRaster <- plotRaster(out_rasterizeGeneExpression_spe)

test_that("Correct structure of rasterizeMatrix output (list)", {
  # output is a list
  expect_type(out_rasterizeMatrix, "list")
  # output has correct element names
  expect_equal(names(out_rasterizeMatrix), c("data_rast", "pos_rast", "meta_rast"))
  # output elements have the same number of pixels
  expect_equal(ncol(out_rasterizeMatrix$data_rast), nrow(out_rasterizeMatrix$pos_rast))
  expect_equal(ncol(out_rasterizeMatrix$data_rast), nrow(out_rasterizeMatrix$meta_rast))
  expect_equal(nrow(out_rasterizeMatrix$pos_rast), nrow(out_rasterizeMatrix$meta_rast))
  # output element (data_rast) has the same number of features (genes or celltypes) as the input
  expect_equal(nrow(out_rasterizeMatrix$data_rast), nrow(data))
  # output element (pos_rast) has x,y coordinates
  expect_equal(colnames(out_rasterizeMatrix$pos_rast), c("x","y"))
  # output number of pixels ≤ number of cells in the input
  expect_lte(ncol(out_rasterizeMatrix$data_rast), ncol(data))
})

test_that("Correct structure of rasterizeGeneExpression output when input is a SpatialExperiment object (SpatialExperiment)", {
  # output is a SpatialExperiment
  expect_s4_class(out_rasterizeGeneExpression_spe, "SpatialExperiment")
  # output has the same number of genes as the input
  
  # output contains an assay slot
  
  # output contains a SpatialCoords slot
  
  # output contains a colData slot
  
})

test_that("Correct structure of rasterizeGeneExpression output when input is a list of SpatialExperiment object (list)", {
  # output is a list
  
  # output list names have the same list names as the input
  
})

test_that("Correct structure of rasterizeCellType output when input is a SpatialExperiment object (SpatialExperiment)", {
  # output is a SpatialExperiment
  
  # output has the same number of celltypes as the input
  
  # output contains an assay slot
  
  # output contains a SpatialCoords slot
  
  # output contains a colData slot
  
})

test_that("Correct structure of rasterizeCellType output when input is a list of SpatialExperiment object (list)", {
  # output is a list
  
  # output list names have the same list names as the input
  
})

test_that("Correct structure of permutateByRotation output (list)", {
  # output is a list
  
  # output has the same number of elements as the n_perm argument
  
  # output does not have redundant rotations
  angles <- as.numeric(gsub("rotated_", "", names(out_permutateByRotation)))
  expect_false(any(duplicated(angles)))
})

test_that("Correct structure of plotRaster output (ggplot)", {
  # output is a ggplot
  expect_s3_class(out_plotRaster, "ggplot")
})