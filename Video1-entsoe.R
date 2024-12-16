
# Tomo ideas https://gist.github.com/Preetam/1f8951294c475eb2188c









library('glasso') 
library('igraph') 
library('qgraph') 
library('MASS') 
library(ggplot2)
library(reshape2)

library(gplots) 

set.seed(100)


#entso = read.csv( '/home/reinaldo/7a310714-2a6d-44bd-bd76-c6a65540eb82/Data Science/Maestria/notebook/Maestria/EEA/Presentacion final/Dataset_Precios/DA.csv', header = TRUE,)
#entso <- as.data.frame(read.csv('/home/reinaldo/7a310714-2a6d-44bd-bd76-c6a65540eb82/Data Science/Maestria/notebook/Maestria/EEA/Presentacion final/Dataset_Precios/DA.csv', header = TRUE))
entso <- read.csv(
  '/home/reinaldo/7a310714-2a6d-44bd-bd76-c6a65540eb82/Data Science/Maestria/notebook/Maestria/EEA/Presentacion final/Dataset_Precios/DA.csv',
  header = TRUE,
  row.names = 1
)


entso <- read.csv(
  '/home/reinaldo/Downloads/DA2.csv',
  header = TRUE,
  row.names = 1
)

########## 8M
start_time <- as.POSIXct("2023-12-31 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
num_rows <- nrow(entso)
timestamps <- seq(from = start_time, by = "hour", length.out = num_rows)
row.names(entso) <- timestamps
entso10= entso[1:10]
entso10$Time <- row.names(entso10)
entso_long <- melt(entso10, id.vars = "Time")


entso_long$Time <- as.POSIXct(entso_long$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
ggplot(entso_long, aes(x = Time, y = value, color = variable)) +
  geom_line() +
  labs(title = "DA Price 2024", x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "right")


############# 1M
start_time <- as.POSIXct("2023-12-31 23:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
num_rows <- nrow(entso)
timestamps <- seq(from = start_time, by = "hour", length.out = num_rows)
row.names(entso) <- timestamps
entso10= entso[1:24*30,1:15]
entso10$Time <- row.names(entso10)
entso_long <- melt(entso10, id.vars = "Time")


entso_long$Time <- as.POSIXct(entso_long$Time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
ggplot(entso_long, aes(x = Time, y = value, color = variable)) +
  geom_line() +
  labs(title = "DA Price 2024", x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "right")




entso_cov = cov(entso)
heatmap.2(
  entso_cov,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  trace = "none",
  col = colorRampPalette(c("white", "red"))(100),
  main = "Heat map Covarianza",
  labRow = colnames(entso),
  labCol = colnames(entso),
  srtCol = 45,  
  cexCol = 0.5,  
  cexRow = 0.5,
  margins = c(8, 8),  
  key = FALSE   
)




adyacencia <- entso_cov> 1E-1; diag(adyacencia) <- 0
adyacencia <- graph_from_adjacency_matrix(adyacencia, mode = 'undirected')
communities <- cluster_fast_greedy(adyacencia)

membership <- communities$membership
order <- order(membership)
sorted_matrix <- as.matrix(adyacencia[order, order])
sorted_row_col_names <- colnames(entso)[order]

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

entso_var<- cov(entso)

g<-glasso(entso_var, rho=0.001)

adyacencia <- abs(g$wi) > 1E-1; diag(adyacencia) <- 0
adyacencia <- graph_from_adjacency_matrix(adyacencia, mode = 'undirected')

communities <- cluster_fast_greedy(adyacencia)

membership <- communities$membership
order <- order(membership)
sorted_matrix <- as.matrix(adyacencia[order, order])
sorted_row_col_names <- colnames(entso)[order]

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
subset_matrix <- sorted_matrix[1:11, 1:11]
subset_row_col_names <- sorted_row_col_names[1:11]

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
  srtCol = 45,  
  cexCol = 1.0,  
  cexRow = 1.0,
  margins = c(12, 12),  
  key = FALSE   
)


#######################

entso_corr = cor(entso)

entso_corr_inv <- solve(entso_corr)
det(entso_corr)



