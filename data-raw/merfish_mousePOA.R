## Prepares data/merfish_mousePOA.rda
## format the MERFISH mouse POA dataset into a SpatialExperiment object for the package

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(here)

data <- read.csv('~/Downloads/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv')

animal <- 1
sex <- "Female"
behavior <- "Naive"
bregma <- "-0.29"
data_sub <- data[(data$Animal_ID == animal & data$Animal_sex == sex & data$Behavior == behavior & data$Bregma == bregma),]
dim(data_sub)

## genes-by-cells matrix
# extract the genes-by-cells matrix as a sparse matrix (dgCMatrix)
mat <- as(t(data_sub[,10:ncol(data_sub)]), "CsparseMatrix")

# remove blank genes used for quality control
blanks <- rownames(mat)[grepl("Blank", rownames(mat))]
mat <- mat[setdiff(rownames(mat),blanks),]

## spatial coordinates matrix
# extract the spatial coordinates
pos <- data_sub[,c("Centroid_X", "Centroid_Y")]
colnames(pos) <- c("x","y")

# make x,y coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## cell-type labels
# extract the data frame with cell-type labels
meta <- data_sub[,c("Bregma", "Cell_class", "Neuron_cluster_ID")]
colnames(meta) <- c("bregma", "celltype", "neurontype")

## standardize cell IDs for the extracted objects
colnames(mat) <- rownames(pos) <- rownames(meta) <- data_sub$Cell_ID

## filter genes with NaN values
bad_genes <- names(which(rowSums(is.nan(mat)) > 0)) 
mat <- mat[setdiff(rownames(mat),bad_genes),]

## filter cells with NaN values
bad_cells <- names(which(colSums(is.nan(mat)) > 0))
mat <- mat[,setdiff(colnames(mat),bad_cells)]
pos <- pos[setdiff(rownames(pos),bad_cells),]
meta <- meta[setdiff(rownames(pos),bad_cells),]

df_plt <- data.frame(pos, total_gexp = colSums(mat))

ggplot(df_plt, aes(x = x, y = y, color = total_gexp)) +
  coord_fixed() +
  geom_point(size = 1, stroke = 0) +
  scale_color_viridis_c(name = "total gene expression") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

df_plt <- data.frame(pos, celltype = meta$celltype)

ggplot(df_plt, aes(x = x, y = y, color = celltype)) +
  coord_fixed() +
  geom_point(size = 1, stroke = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

merfish_mousePOA <- SpatialExperiment::SpatialExperiment(
  assays = list(volnorm = mat),
  spatialCoords = as.matrix(pos),
  colData = meta
)

usethis::use_data(merfish_mousePOA, overwrite = TRUE)
