library(flowCore)
library(FlowSOM)
library(ConsensusClusterPlus)
library(RColorBrewer)
library(pheatmap)
library(dplyr)

setwd("/Users/monikagrabowska/Downloads/FlowRepository_FR-FCM-ZZKZ_files/")

fcs <- read.FCS("AML_normal_viSNEgates_concat.fcs", transformation = FALSE, truncate_max_range = FALSE)
param <- pData(parameters(fcs))
param$desc <- gsub("-1.*", "", param$desc)
param$desc <- gsub("-", "", param$desc)

lineage_markers <- param$desc[5:31]

colnames(fcs) <- param$desc

fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

mst <- BuildMST(som)

metaClustering1 <- as.character(metaClustering_consensus(mst$map$codes,k = 3))
metaClustering2 <- as.character(metaClustering_consensus(mst$map$codes,k = 12))
metaClustering3 <- as.character(metaClustering_consensus(mst$map$codes,k = 90))

PlotStars(mst, backgroundValues = metaClustering1)
PlotStars(mst, backgroundValues = metaClustering2)
PlotStars(mst, backgroundValues = metaClustering3)

# query to look for CD45+ cells
query <- c("CD45" = "high")
query_res <- QueryStarPlot(mst, query, equalNodeSize = TRUE, plot = FALSE)
cellTypes <- factor(rep("Unlabeled", mst$map$nNodes), levels=c("Unlabeled", "CD45 cells"))
cellTypes[query_res$selected] <- "CD45 cells"
p <- PlotStars(mst, backgroundValues = cellTypes, backgroundColor = c("#FFFFFF00","#0000FF"))

# metaclustering into 20 clusters 
codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)

# get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

plot_clustering_heatmap_wrapper <- function(expr, cell_clustering, color_clusters, cluster_merging = NULL) {
  
  # calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  
  # calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr_median[, colnames(expr)])
  rownames(expr_heat) <- expr_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  breaksList = c(seq(0, 9, by = 1),
                 seq(10, 99, by = 10),
                 seq(100, 1000, by = 100))
  
  # colors for the heatmap
  color <- colorRampPalette(c("black", "yellow"))(length(breaksList))
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = FALSE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = FALSE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend, breaks = breaksList,
           legend_breaks = c(0,1,10,100,1000), legend_labels = c("0","","","100","1000"))
}

plot_clustering_heatmap_wrapper(expr = exprs(fcs[,lineage_markers]), cell_clustering = cell_clustering1, color_clusters = color_clusters)
