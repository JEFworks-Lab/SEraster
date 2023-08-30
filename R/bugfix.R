library(SEraster)
library(Matrix)
library(raster)
library(ggplot2)

## Jean's hard coding for now (replace with demo dataset to be included with package)
meta <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Pallav_mouseheart_MERFISH/Sham1-Newseg/Sham1_Newseg_afterQC_cell_meta.csv', row.names=1)
gexp <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Pallav_mouseheart_MERFISH/Sham1-Newseg/Sham1_Newseg_afterQC_cell_by_gene.csv', row.names=1)
pos <- meta[, c('center_x', 'center_y')]
mat <- t(Matrix(as.matrix(gexp/meta$volume), sparse=TRUE))
class(mat)
dim(pos)
dim(mat)

## this produces an error
# out <- SEraster::rasterizeSparseMatrix(mat, pos, fun="mean")

resolution = 100
fun = "mean"
data = mat
## create RasterLayer (simplified way, using actual operations in raster() function)

ext <- raster::extent(min(pos[,1]), max(pos[,1]), min(pos[,2]), max(pos[,2]))
r <- methods::new('RasterLayer', extent=ext)
raster::res(r) <- resolution
r

ncols <- round((max(pos[,1]) - min(pos[,1]))/resolution)
nrows <- round((max(pos[,2]) - min(pos[,2]))/resolution)
r <- raster::raster(ncols=ncols, nrows=nrows,
                    xmn=min(pos[,1]), ymn=min(pos[,2]), xmx=max(pos[,1]), ymx=max(pos[,2]))
## if you don't set resolution, raster::raster adjusts the resolution to fit all points
## probably better to set exact resolution and add padding to prevent some points from being set out of the bounding box
r

## convert RasterLayer to SpatialPoints (sp package)
## might be good to use sf package since sp package is retiring soon
pts <- raster::rasterToPoints(r, spatial = TRUE)
## get x,y coordinates for each pixel
pos_pixel <- pts@coords
rownames(pos_pixel) <- paste0("pixel", 1:nrow(pos_pixel))

## assign pixel ID to each single cell xy coordinates
pixel_ids <- raster::cellFromXY(r, pos)
length(pixel_ids)
names(pixel_ids) <- rownames(pos)

## NAs in these pixel_ids
table(is.na(pixel_ids))

## plot cells without any pixel IDs
pixel_ids[is.na(pixel_ids)]

cell_ids <- names(pixel_ids[is.na(pixel_ids)])

df <- data.frame(pos)
df2 <- df[cell_ids,]

ggplot(df, aes(x = center_x, y = center_y)) +
  geom_point(size = 0.1, color = "lightgray") +
  geom_point(data = df2, aes(x = center_x, y = center_y), size = 0.1, color = "red") +
  theme_classic()

r
