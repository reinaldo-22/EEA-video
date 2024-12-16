

library('glasso') 
library('igraph') 
library('qgraph') 
library('MASS') 
library(ggplot2)
library(reshape2)

library(gplots) 

set.seed(100)



library(abind)
array_2d <- read.csv('/home/reinaldo/Downloads/sq_small3.csv', header = FALSE)


array_3d <- abind(array_2d, along = 3)
dim(array_3d)


#dim(array_3d) <- c(1590, 350,160)
#dim(array_3d) <- c(901, 101, 61)
#dim(array_3d) <- c(901, 61, 61)
dim(array_3d) <- c(901, 61, 61)


lat_names <- paste0("lat#", 1:dim(array_3d)[2])
lon_names <- paste0("lon#", 1:dim(array_3d)[3])

dimnames(array_3d) <- list(NULL, lat_names, lon_names)

slice_1 <- array_3d[50,,]
dim(slice_1)
slice_1 <- slice_1[, ncol(slice_1):1] 

image(slice_1, main = "Slice at Axis = 1", xlab = "Dimension 2", ylab = "Dimension 3", col = heat.colors(256))
#subsampled_array= array_3d
subsampled_array <- array_3d[, seq(1, dim(array_3d)[2], by = 3), seq(1, dim(array_3d)[3], by = 3)]
dimnames(subsampled_array) <- list(
  NULL,
  lat_names[seq(1, length(lat_names), by = 3)],
  lon_names[seq(1, length(lon_names), by = 3)]
)


dim(subsampled_array)

slice_1 <- subsampled_array[50,,]
dim(slice_1)
slice_1 <- slice_1[, ncol(slice_1):1] 

image(slice_1, main = "Slice at Axis = 1", xlab = "Dimension 2", ylab = "Dimension 3", col = heat.colors(256))


#dim(subsampled_array) <- c(1590, dim(subsampled_array)[2]*dim(subsampled_array)[3])
dim(subsampled_array) <- c(901, dim(subsampled_array)[2]*dim(subsampled_array)[3])

dim( subsampled_array)




subsampled_array_cov = cov(subsampled_array)
lat_sub <- lat_names[seq(1, length(lat_names), by = 3)]
lon_sub <- lon_names[seq(1, length(lon_names), by = 3)]
dim_names <- paste0(rep(lat_sub, each = length(lon_sub)), "_", rep(lon_sub, times = length(lat_sub)))

dimnames(subsampled_array_cov) <- list(dim_names, dim_names)

heatmap(
  subsampled_array_cov,
  Rowv = NA,
  Colv = NA,
  col = colorRampPalette(c("white", "red"))(100),
  main = "Heat Map Covariance",
  labRow = NA,
  labCol = NA,
  margins = c(8, 8)
)

heatmap.2(
  subsampled_array_cov,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  trace = "none",
  col = colorRampPalette(c("white", "red"))(100),
  main = "Heat map Covarianza",
  labRow = colnames(subsampled_array_cov),
  labCol = colnames(subsampled_array_cov),
  srtCol = 45,  
  cexCol = 0.5,  
  cexRow = 0.5,
  margins = c(8, 8),  
  key = FALSE   
)




adyacencia <- subsampled_array_cov> 1E-1
diag(adyacencia) <- 0
adyacencia <- graph_from_adjacency_matrix(adyacencia, mode = 'undirected')


communities <- cluster_fast_greedy(adyacencia)

membership <- communities$membership
order <- order(membership)
sorted_matrix <- as.matrix(adyacencia[order, order])
sorted_row_col_names <- colnames(subsampled_array_cov)[order]

heatmap.2(
  sorted_matrix,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  trace = "none",
  col = colorRampPalette(c("white", "red"))(100),
  main = "Heat map Covarianza con comunidades",
  labRow = sorted_row_col_names,
  labCol = sorted_row_col_names,
  srtCol = 45,  
  cexCol = 0.5,  
  cexRow = 0.5,
  margins = c(8, 8),  
  key = FALSE   
)



############## Glasso adyacencia

#entso_var<- cov(entso)


g<-glasso(subsampled_array_cov, rho=0.15, maxit=50)

adyacencia <- abs(g$wi) > 1E-1; diag(adyacencia) <- 0
adyacencia <- graph_from_adjacency_matrix(adyacencia, mode = 'undirected')

communities <- cluster_fast_greedy(adyacencia)

membership <- communities$membership
order <- order(membership)
sorted_matrix <- as.matrix(adyacencia[order, order])
sorted_row_col_names <- colnames(subsampled_array_cov)[order]

heatmap.2(
  sorted_matrix,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  trace = "none",
  col = colorRampPalette(c("white", "red"))(100),
  main = "Grafo con comunidades",
  labRow = sorted_row_col_names,
  labCol = sorted_row_col_names,
  srtCol = 45,  
  cexCol = 0.5,  
  cexRow = 0.5,
  margins = c(8, 8),  
  key = FALSE   
)


####### Sub seleccion de comunidad
subset_matrix <- sorted_matrix[1:20, 1:20]
subset_row_col_names <- sorted_row_col_names[1:20]

heatmap.2(
  subset_matrix,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  trace = "none",
  col = colorRampPalette(c("white", "red"))(100),
  main = "Grafo con comunidades",
  labRow = subset_row_col_names,
  labCol = subset_row_col_names,
  srtCol = 90,  
  cexCol = 0.7,  
  cexRow = 0.7,
  margins = c(12, 12),  
  key = FALSE   
)


#######################

det(subsampled_array_cov)


