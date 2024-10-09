# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# 
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'))
# 
# install.packages("devtools")
# devtools::install_github('cole-trapnell-lab/monocle3')


# install.packages("Seurat")
# devtools::install_github("satijalab/seurat-wrappers")
# remotes::install_github("mojaveazure/seurat-disk")

# Libraries --------------------------------------------------------------------

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(zellkonverter)
library(SeuratWrappers)
library(data.table)
library(tidyverse)
library(scales)

# Data loading and creating cds object -----------------------------------------
sce1=readH5AD("../adata_macrophages.h5ad", verbose = TRUE)
adata_Seurat <- as.Seurat(sce1, counts = "X"
                          , data = NULL
)

rm(sc)

cds <- as.cell_data_set(adata_Seurat)

rm(adata_Seurat)

# Processing data --------------------------------------------------------------

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = F)

## Step 6: Order cells
cds <- order_cells(cds)

# Data visualisation -----------------------------------------------------------

# Single cell trajectory plot
plot_cells(cds, show_trajectory_graph = F)

# Pseudotime single cell trajectory plot
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# Dataframe wrangling for barplots ---------------------------------------------

clusters <- as.data.frame(cds@clusters@listData$UMAP$clusters)

cds@colData$sample_new <- sapply(cds@colData$Sample, function(row) {
  if (row == "GSM6593315") {
    return("AD1")
  } else if (row == "GSM6593316") {
    return("AD2")
  } else if (row == "GSM6593317") {
    return("AD3")
  } else if (row == "GSM6593318") {
    return("AD4")
  } else if (row == "GSM6593319") {
    return("AD5")
  } else if (row == "GSM6593320") {
    return("AD6")
  } else if (row == "GSM6593321") {
    return("Normal1")
  } else if (row == "GSM6593322") {
    return("Normal2")
  } else if (row == "GSM6593323") {
    return("Normal3")
  } else {
    return(NA) # Add an NA for unmatched cases
  }
})

sample_new <- as.data.frame(cds@colData)[c("cell_id", "sample_new")]

clusters$cell_id <- row.names(clusters)

new_df <- merge(clusters, sample_new, by="cell_id")
new_df
colnames(new_df) <- c("cell_id", "clusters", "sample_new")

new_df_freq <- new_df %>% group_by(sample_new, clusters) %>% 
  summarise(freq = n())

new_df_freq$condition <- new_df_freq$sample_new %>% sapply(function(row) {
  if (substr(row, start=1, stop=2) == "AD"){
    return("AD")
  } else {
    return("Normal")
  }
})

# Barplotting ==================================================================

# Barplot of the AD and normal samples
ggplot(new_df_freq, aes(fill=clusters, x=sample_new, y=freq)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab(NULL) +
  theme_classic() +
  ylab("Percentage of cells") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(labels = percent, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) 

# Barplot of the patients samples grouped either to AD or normal group
ggplot(new_df_freq, aes(fill=clusters, x=condition, y=freq)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab(NULL) +
  theme_classic() +
  ylab("Percentage of cells") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(labels = percent, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) 
