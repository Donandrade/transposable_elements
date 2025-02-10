# Install needed packages
install.packages("rentrez", dependencies = TRUE)
install.packages("devtools")
install.packages("ape")
install.packages("seqinr")

# Load required libraries
library(seqinr)
library(dendextend)
library(circlize)
library(ape)
library(cluster)


# Read alignment data
shrooms_seq <- read.alignment("mule_MuDR_mafft_trimed/mule_MuDR_mafft_trimed.fasta", format = "fasta")

# Calculate distance matrix
shrooms_dist <- dist.alignment(shrooms_seq, matrix = "identity")

shrooms_dist[is.na(shrooms_dist)] <- mean(shrooms_dist, na.rm = TRUE)

# Perform K-means clustering
clusters <- kmeans(shrooms_dist, 20)

cluster_matrix <- as.matrix(clusters$cluster)


colnames(cluster_matrix)

write.table(cluster_matrix, file = "caulimovirus_mafft_trimed/clusters_genes_caulimovirus.csv", quote = FALSE, sep = "\t")

# Calcular média das distâncias para cada cluster
avg_distances <- numeric(20)

paste("caulimovirus_mafft_trimed/","cluster_", "1", "_caulimovirus_distance.csv", sep = "")

# Inicializar uma matriz vazia para armazenar os resultados
resultados <- matrix(0, nrow = 20, ncol = 5) # 20 linhas para os clusters, 4 colunas para os resultados

for (i in 1:20) {
  cluster_indices <- which(clusters$cluster == i)
  file_name <- paste("caulimovirus_mafft_trimed/","cluster", i, "_caulimovirus_distance.csv", sep = "")
  write.csv(cluster_indices, file = file_name)
  cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]
  avg_distances[i] <- mean(cluster_distances)
  tmp <- c(rownames(as.table(which(clusters$cluster == i))))
  conta_cc <- 0
  conta_CA <- 0
  conta_CB <- 0
  for (x in 1:length(tmp)){
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CC"){
      conta_cc <- conta_cc+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CA"){
      conta_CA <- conta_CA+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CB"){
      conta_CB <- conta_CB+1
    }
  }
  
  # Preencher a matriz com os resultados
  resultados[i, 1] <- i
  resultados[i, 2] <- avg_distances[i]
  resultados[i, 3] <- conta_cc
  resultados[i, 4] <- conta_CA
  resultados[i, 5] <- conta_CB
}

# Definir nomes para as colunas
colnames(resultados) <- c("Cluster number", "Avg_Distances", "CC sequences number", 
                          "CA sequences number", "CB sequences number")

# Imprimir a matriz resultante
write.table(resultados, file = "caulimovirus_mafft_trimed/media_distancias_contagemSeq_caulimovirus.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Find cluster with the minimum average distance
cluster_with_min_distance <- which.min(avg_distances)

cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Find cluster with the seccond minimum average distance
nova_lista <- subset(avg_distances, avg_distances != avg_distances[cluster_with_min_distance])

cluster_with_min_distance <- which.min(nova_lista)
nova_lista
cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Generate dendrogram
shrooms_hclust <- hclust(shrooms_dist, method = "complete")
hc <- as.dendrogram(shrooms_hclust)
phy <- as.phylo(hc)

save.image(file = "caulimovirus_mafft_trimed/caulimovirus.RData")
write.tree(phy, file = "caulimovirus_mafft_trimed/caulimovirus_dendrograma.nwk")

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Read alignment data
shrooms_seq <- read.alignment("line_L1_mafft_trimed/line_L1_mafft_trimed.fasta", format = "fasta")

# Calculate distance matrix
shrooms_dist <- dist.alignment(shrooms_seq, matrix = "identity")

shrooms_dist[is.na(shrooms_dist)] <- mean(shrooms_dist, na.rm = TRUE)

# Perform K-means clustering
clusters <- kmeans(shrooms_dist, 20)

cluster_matrix <- as.matrix(clusters$cluster)

write.table(cluster_matrix, file = "line_L1_mafft_trimed/clusters_genes_line_L1.csv", quote = FALSE, sep = "\t")

# Calcular média das distâncias para cada cluster
avg_distances <- numeric(20)

# Inicializar uma matriz vazia para armazenar os resultados
resultados <- matrix(0, nrow = 20, ncol = 5) # 20 linhas para os clusters, 4 colunas para os resultados

for (i in 1:20) {
  cluster_indices <- which(clusters$cluster == i)
  file_name <- paste("line_L1_mafft_trimed/","cluster", i, "_line_L1_distance.csv", sep = "")
  write.csv(cluster_indices, file = file_name)
  cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]
  avg_distances[i] <- mean(cluster_distances)
  tmp <- c(rownames(as.table(which(clusters$cluster == i))))
  conta_cc <- 0
  conta_CA <- 0
  conta_CB <- 0
  for (x in 1:length(tmp)){
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CC"){
      conta_cc <- conta_cc+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CA"){
      conta_CA <- conta_CA+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CB"){
      conta_CB <- conta_CB+1
    }
  }
  
  # Preencher a matriz com os resultados
  resultados[i, 1] <- i
  resultados[i, 2] <- avg_distances[i]
  resultados[i, 3] <- conta_cc
  resultados[i, 4] <- conta_CA
  resultados[i, 5] <- conta_CB
}

# Definir nomes para as colunas
colnames(resultados) <- c("Cluster number", "Avg_Distances", "CC sequences number", 
                          "CA sequences number", "CB sequences number")

# Imprimir a matriz resultante
write.table(resultados, file = "line_L1_mafft_trimed/media_distancias_contagemSeq_line_L1.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Find cluster with the minimum average distance
cluster_with_min_distance <- which.min(avg_distances)

cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Find cluster with the seccond minimum average distance
nova_lista <- subset(avg_distances, avg_distances != avg_distances[cluster_with_min_distance])

cluster_with_min_distance <- which.min(nova_lista)
nova_lista
cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Generate dendrogram
shrooms_hclust <- hclust(shrooms_dist, method = "complete")
hc <- as.dendrogram(shrooms_hclust)
phy <- as.phylo(hc)

save.image(file = "line_L1_mafft_trimed/line_L1.RData")
write.tree(phy, file = "line_L1_mafft_trimed/line_L1_dendrograma.nwk")

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Read alignment data
shrooms_seq <- read.alignment("line_RTE_mafft_trimed/line_RTE_mafft_trimed.fasta", format = "fasta")

# Calculate distance matrix
shrooms_dist <- dist.alignment(shrooms_seq, matrix = "identity")

shrooms_dist[is.na(shrooms_dist)] <- mean(shrooms_dist, na.rm = TRUE)

# Perform K-means clustering
clusters <- kmeans(shrooms_dist, 20)

cluster_matrix <- as.matrix(clusters$cluster)

write.table(cluster_matrix, file = "line_RTE_mafft_trimed/clusters_genes_line_RTE.csv", quote = FALSE, sep = "\t")

# Calcular média das distâncias para cada cluster
avg_distances <- numeric(20)

# Inicializar uma matriz vazia para armazenar os resultados
resultados <- matrix(0, nrow = 20, ncol = 5) # 20 linhas para os clusters, 4 colunas para os resultados

for (i in 1:20) {
  cluster_indices <- which(clusters$cluster == i)
  file_name <- paste("line_RTE_mafft_trimed/","cluster", i, "_line_RTE_distance.csv", sep = "")
  write.csv(cluster_indices, file = file_name)
  cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]
  avg_distances[i] <- mean(cluster_distances)
  tmp <- c(rownames(as.table(which(clusters$cluster == i))))
  conta_cc <- 0
  conta_CA <- 0
  conta_CB <- 0
  for (x in 1:length(tmp)){
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CC"){
      conta_cc <- conta_cc+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CA"){
      conta_CA <- conta_CA+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CB"){
      conta_CB <- conta_CB+1
    }
  }
  
  # Preencher a matriz com os resultados
  resultados[i, 1] <- i
  resultados[i, 2] <- avg_distances[i]
  resultados[i, 3] <- conta_cc
  resultados[i, 4] <- conta_CA
  resultados[i, 5] <- conta_CB
}

# Definir nomes para as colunas
colnames(resultados) <- c("Cluster number", "Avg_Distances", "CC sequences number", 
                          "CA sequences number", "CB sequences number")

# Imprimir a matriz resultante
write.table(resultados, file = "line_RTE_mafft_trimed/media_distancias_contagemSeq_line_RTE.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Find cluster with the minimum average distance
cluster_with_min_distance <- which.min(avg_distances)

cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Find cluster with the seccond minimum average distance
nova_lista <- subset(avg_distances, avg_distances != avg_distances[cluster_with_min_distance])

cluster_with_min_distance <- which.min(nova_lista)
nova_lista
cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Generate dendrogram
shrooms_hclust <- hclust(shrooms_dist, method = "complete")
hc <- as.dendrogram(shrooms_hclust)
phy <- as.phylo(hc)

save.image(file = "line_RTE_mafft_trimed/line_RTE.RData")
write.tree(phy, file = "line_RTE_mafft_trimed/line_RTE_dendrograma.nwk")

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Read alignment data
shrooms_seq <- read.alignment("mule_MuDR_mafft_trimed/mule_MuDR_mafft_trimed.fasta", format = "fasta")

# Calculate distance matrix
shrooms_dist <- dist.alignment(shrooms_seq, matrix = "identity")

shrooms_dist[is.na(shrooms_dist)] <- mean(shrooms_dist, na.rm = TRUE)

# Perform K-means clustering
clusters <- kmeans(shrooms_dist, 20)

cluster_matrix <- as.matrix(clusters$cluster)

write.table(cluster_matrix, file = "mule_MuDR_mafft_trimed/clusters_genes_mule_MuDR.csv", quote = FALSE, sep = "\t")

# Calcular média das distâncias para cada cluster
avg_distances <- numeric(20)

# Inicializar uma matriz vazia para armazenar os resultados
resultados <- matrix(0, nrow = 20, ncol = 5) # 20 linhas para os clusters, 4 colunas para os resultados

for (i in 1:20) {
  cluster_indices <- which(clusters$cluster == i)
  file_name <- paste("mule_MuDR_mafft_trimed/","cluster", i, "_mule_MuDR_distance.csv", sep = "")
  write.csv(cluster_indices, file = file_name)
  cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]
  avg_distances[i] <- mean(cluster_distances)
  tmp <- c(rownames(as.table(which(clusters$cluster == i))))
  conta_cc <- 0
  conta_CA <- 0
  conta_CB <- 0
  for (x in 1:length(tmp)){
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CC"){
      conta_cc <- conta_cc+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CA"){
      conta_CA <- conta_CA+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CB"){
      conta_CB <- conta_CB+1
    }
  }
  
  # Preencher a matriz com os resultados
  resultados[i, 1] <- i
  resultados[i, 2] <- avg_distances[i]
  resultados[i, 3] <- conta_cc
  resultados[i, 4] <- conta_CA
  resultados[i, 5] <- conta_CB
}

# Definir nomes para as colunas
colnames(resultados) <- c("Cluster number", "Avg_Distances", "CC sequences number", 
                          "CA sequences number", "CB sequences number")

# Imprimir a matriz resultante
write.table(resultados, file = "mule_MuDR_mafft_trimed/media_distancias_contagemSeq_mule_MuDR.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Find cluster with the minimum average distance
cluster_with_min_distance <- which.min(avg_distances)

cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Find cluster with the seccond minimum average distance
nova_lista <- subset(avg_distances, avg_distances != avg_distances[cluster_with_min_distance])

cluster_with_min_distance <- which.min(nova_lista)
nova_lista
cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Generate dendrogram
shrooms_hclust <- hclust(shrooms_dist, method = "complete")
hc <- as.dendrogram(shrooms_hclust)
phy <- as.phylo(hc)

save.image(file = "mule_MuDR_mafft_trimed/mule_MuDR.RData")
write.tree(phy, file = "mule_MuDR_mafft_trimed/mule_MuDR_dendrograma.nwk")

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Read alignment data
shrooms_seq <- read.alignment("sine_mafft_trimed/sine_mafft_trimed.fasta", format = "fasta")

# Calculate distance matrix
shrooms_dist <- dist.alignment(shrooms_seq, matrix = "identity")

shrooms_dist[is.na(shrooms_dist)] <- mean(shrooms_dist, na.rm = TRUE)

# Perform K-means clustering
clusters <- kmeans(shrooms_dist, 20)

cluster_matrix <- as.matrix(clusters$cluster)

write.table(cluster_matrix, file = "sine_mafft_trimed/clusters_genes_sine.csv", quote = FALSE, sep = "\t")

# Calcular média das distâncias para cada cluster
avg_distances <- numeric(20)

# Inicializar uma matriz vazia para armazenar os resultados
resultados <- matrix(0, nrow = 20, ncol = 5) # 20 linhas para os clusters, 5 colunas para os resultados

for (i in 1:20) {
  cluster_indices <- which(clusters$cluster == i)
  file_name <- paste("sine_mafft_trimed/","cluster", i, "_sine_distance.csv", sep = "")
  write.csv(cluster_indices, file = file_name)
  cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]
  avg_distances[i] <- mean(cluster_distances)
  tmp <- c(rownames(as.table(which(clusters$cluster == i))))
  conta_cc <- 0
  conta_CA <- 0
  conta_CB <- 0
  for (x in 1:length(tmp)){
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CC"){
      conta_cc <- conta_cc+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CA"){
      conta_CA <- conta_CA+1
    }
    if(substr(tmp[x], nchar(tmp[x]) - 1, nchar(tmp[x])) == "CB"){
      conta_CB <- conta_CB+1
    }
  }
  
  # Preencher a matriz com os resultados
  resultados[i, 1] <- i
  resultados[i, 2] <- avg_distances[i]
  resultados[i, 3] <- conta_cc
  resultados[i, 4] <- conta_CA
  resultados[i, 5] <- conta_CB
}

# Definir nomes para as colunas
colnames(resultados) <- c("Cluster number",
                          "Avg_Distances",
                          "CC sequences number", 
                          "CA sequences number",
                          "CB sequences number")

# Imprimir a matriz resultante
write.table(resultados, file = "sine_mafft_trimed/media_distancias_contagemSeq_sine.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Find cluster with the minimum average distance
cluster_with_min_distance <- which.min(avg_distances)

cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Find cluster with the seccond minimum average distance
nova_lista <- subset(avg_distances, avg_distances != avg_distances[cluster_with_min_distance])

cluster_with_min_distance <- which.min(nova_lista)
nova_lista
cluster_with_min_distance

cluster_indices <- which(clusters$cluster == cluster_with_min_distance)

cluster_distances <- as.matrix(shrooms_dist)[cluster_indices, cluster_indices]

cluster_hclust <- hclust(as.dist(cluster_distances), method = "complete")

hc <- as.dendrogram(cluster_hclust)

circlize_dendrogram(hc, dend_track_height = 0.6, labels = TRUE, labels_size = 20)


# Generate dendrogram
shrooms_hclust <- hclust(shrooms_dist, method = "complete")
hc <- as.dendrogram(shrooms_hclust)
phy <- as.phylo(hc)

save.image(file = "sine_mafft_trimed/sine.RData")
write.tree(phy, file = "sine_mafft_trimed/sine_dendrograma.nwk")

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

