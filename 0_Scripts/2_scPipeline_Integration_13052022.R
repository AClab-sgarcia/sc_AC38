## Date: 13/05/2022

## In this pipeline this steps are done: 
##  1. Pre-processing workflow - Normalization + FindVariableFeatures + ScaleData OR SCT 
##  2. Integration
##  3. Linear dimensional reduction
##  4. Clustering
##  5. Marker identification

## Directory structure:
##  1_SeuratObject 
##  2_QC_Figures
##  3_Integration_Figures

dir.create("3_Integration_Figures")

## General pipeline for single-cell data analysis
## Following:
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Guided%20Clustering%20Tutorial%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Introduction%20to%20scRNA-seq%20integration%20(May2022).R 
## -  Single-cell RNA-seq data analysis workshop: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html 

## Definitions of the functions from: https://cloud.r-project.org/web/packages/Seurat/Seurat.pdf (Version 4.1.1 || Date 2022-05-01)

####################################################################
## Libraries
####################################################################
## librearies that need to be loaded in nextera docker
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("dplyr"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("purrr"))
suppressPackageStartupMessages(require("future"))
source("utils/custom_seurat_functions.R")

# suppressPackageStartupMessages(require("scater"))
# suppressPackageStartupMessages(require("patchwork"))
# suppressPackageStartupMessages(require("cowplot"))
# suppressPackageStartupMessages(require("celldex"))
# suppressPackageStartupMessages(require("SingleR"))
# suppressPackageStartupMessages(require("scales"))
# suppressPackageStartupMessages(require("SeuratWrappers"))
# # remotes::install_github("satijalab/seurat-wrappers")
# suppressPackageStartupMessages(require("biomaRt"))
# suppressPackageStartupMessages(require("SeuratDisk"))
# # remotes::install_github("mojaveazure/seurat-disk")

###################################################################
## Variable prep
###################################################################
# The variables are stored in the log file from the 1_scPipeline_QC
log_file <- read.delim("1_scPipeline_QC_13052022_LOG.log", header = TRUE)

project_name <- log_file$Project_name
samples_v <- log_file$Samples
organism <- log_file$Organism
approach <- log_file$Approach
filterting <- log_file$Filtering_parameters
pc_number_doublet <- log_file$PC_number_doublet

print(paste("According to your log file, you are currently usign: ", approach, sep = ""))
print("If your approach is Norm_Feature_Scale, please choose rPCA or CCA for integration. CCA it's more accurate when cell types are conserved but could lead to overcorrection when a large proportion of cells are non overlaping across datasets and it is more computationally expensive. (https://satijalab.org/seurat/articles/integration_rpca.html")

samples_v <- unlist(strsplit(samples_v,","))

# Integration type choose CCA or PCA for the integration (options: "rPCA" || "CCA")
print("Remember to choose between rPCA or CCA for the integration.")
integration_reduction <- "CCA"

for (sample in samples_v){
  # We are going to save each of the objects separately 
  RDS <- readRDS(paste("1_SeuratObjects/QC_", sample, "_", project_name, ".rds", sep = ""))
  column_name <- colnames(RDS@meta.data)[grep("pANN", colnames(RDS@meta.data))]
  RDS[[column_name]] <- NULL
  assign(sample, RDS)
}

if (length(samples_v) > 1){
  data_combined <- merge(x = get(samples_v[1]), y = sapply(samples_v[-1], get), add.cell.ids = samples_v, project = project_name)
}

###################################################################
## 5. Pre-processing workflow
###################################################################
# There are 2 approaches: 
# 1. Normalization + FindVariableFeatures + ScaleData workflow
# 2. SCT workflow

# Iterating over samples in a dataset
# Since we have several samples in our dataset (from different conditions), we want to keep them as separate objects and transform them as that is what is required for integration.
data_combined <- SetIdent(data_combined, value = "orig.ident")

split_data <- SplitObject(data_combined, split.by = "orig.ident")

options(future.globals.maxSize = 2000 * 1024^2)
plan("multicore", workers = 12)

if (approach == "Norm_Feature_Scale"){
  
  # https://satijalab.org/seurat/articles/integration_introduction.html
  
  # normalize and identify variable features for each dataset independently
  split_data <- lapply(X = split_data, FUN = function(x) {
    ##################################################
    ## Normalizing the data
    ##################################################
    # By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in data_combined[["RNA"]]@data.
    # normalization.method: Method for normalization. There are three: LogNormalize, CLR, RC. Default: normalization.method = "LogNormalize"
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    
    ##################################################
    ## Identification of highly variable features (feature selection)
    ##################################################
    # We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
    
    # FindVariableFeatures: identifies features that are outliers on a ’mean variability plot’.
    # Selection.method: How to choose top variable features. There are three: vst, mean.var.plot, dispersion. Default: vst
    # nfeatures = 2000 by default
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # Identify the 15 most highly variable genes
  top15 <- head(VariableFeatures(split_data), 15)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(split_data)
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
  plot2
  
} else if (approach == "SCT"){
  ##################################################
  ## SCTransform - https://satijalab.org/seurat/articles/sctransform_vignette.html
  ##################################################
  
  # Now we will use a ‘for loop’ to run the SCTransform() on each sample, and regress out mitochondrial expression by specifying in the vars.to.regress argument of the SCTransform() function.
  # Before we run this for loop, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:
  options(future.globals.maxSize = 4000 * 1024^3)
  
  # SCTransform: This function calls sctransform::vst. The sctransform package is available at https://github.com/ChristophH/sctransform. Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData work-flow. Results are saved in a new assay (named SCT by default) with counts being (corrected) counts, data being log1p(counts), scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay.
  # Now, we run the following loop to perform the sctransform on all samples. This may take some time (~10 minutes):
  split_data <- lapply(X = split_data, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("percent_mt"))
    })
}  

###################################################################
## 2. Integration
###################################################################
# First, we need to specify that we want to use all of the 3000 most variable genes identified by SCTransform for the integration. By default, this function only selects the top 2000 genes.
# SelectIntegrationFeatures: Choose the features to use when integrating multiple datasets. This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. It returns the top scoring features by this ranking. 
# https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
integ_features <- SelectIntegrationFeatures(object.list = split_data, nfeatures = 3000) 

if (approach == "Norm_Feature_Scale"){
  
  if(integration_reduction == "CCA"){
    integ_anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = integ_features, reduction = "cca", dim = 1:50)
    data_integrated <- IntegrateData(anchorset = integ_anchors, dim = 1:50)
  }
  
  if(integration_reduction == "rPCA"){
    split_data <- lapply(X = split_data, FUN = function(x) {
      ##################################################
      ## Scaling the data
      ##################################################
      # We apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
      # The ScaleData() function: 
      #   1. Shifts the expression of each gene, so that the mean expression across cells is 0
      #   2. Scales the expression of each gene, so that the variance across cells is 1 (gives equal weight in downstream analyses, so that highly-expressed 
      #       genes do not dominate)
      #   3. The results of this are stored in s_01_data[["RNA"]]@scale.data
      x <- ScaleData(x, features = integ_features, verbose = FALSE)
      x <- RunPCA(x, features = integ_features, verbose = FALSE)
    })
    
    # Perform CCA, find the best anchors and filter incorrect anchors. (Note: the progress bar in your console will stay at 0%, but know that it is actually running.)
    # FindIntegrationAnchors: Find a set of anchors between a list of Seurat objects. 
    # normalization.method: LogNormalize or SCT
    integ_anchors <- FindIntegrationAnchors(object.list = split_data, normalization.method = "LogNormalize", anchor.features = integ_features, reduction = "rpca")
    
    # Integrate across conditions.
    # IntegrateData: Perform dataset integration using a pre-computed AnchorSet.
    data_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "LogNormalize")
  }
  
  data_integrated <- ScaleData(data_integrated, verbose = FALSE)
    
} else if (approach == "SCT"){
  # Now, we need to prepare the SCTransform object for integration.
  # PrepSCTIntegration: Prepare an object list normalized with sctransform for integration.
  # Prepare the SCT list object for integration
  split_data <- PrepSCTIntegration(object.list = split_data, anchor.features = integ_features)
  split_data <- lapply(X = split_data, FUN = RunPCA, features = integ_features)
  integ_anchors <- FindIntegrationAnchors(object.list = split_data, normalization.method = "SCT", anchor.features = integ_features, reduction = tolower(integration_reduction))
  data_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
}

####################################################################
## 3. Linear dimensional reduction
####################################################################
DefaultAssay(data_integrated) <- "integrated"

# Run PCA on the scaled data. By default, only the previously determined variable features are used as input, it is possible to chose a certain subset.
# RunPCA: Run a PCA dimensionality reduction. 
data_integrated <- RunPCA(object = data_integrated)

# Plot PCA
PCAPlot(data_integrated, group.by = "orig.ident", split.by = "orig.ident")  
ggsave(paste("3_Integration_Figures/PCAPlot_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

# From https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- data_integrated[["pca"]]@stdev / sum(data_integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
max_dim <- min(co1, co2)

# Run UMAP
data_integrated <- RunUMAP(data_integrated, dims = 1:max_dim, reduction = "pca", metric = "euclidean", min.dist = 0.6)

# VizDimLoadings: Visualize top genes associated with reduction components
VizDimLoadings(data_integrated, dims = 1:2, reduction = "pca")
ggsave(paste("3_Integration_Figures/VizDimLoadings_top_genes_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Plot UMAP
# DimPlot: Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a cell and it’s positioned based on the cell embeddings determined by the reduction technique. 
DimPlot(data_integrated, group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/Integration_UMAP_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Plot UMAP split by sample (Side-by-side comparison of clusters)
# Sometimes it’s easier to see whether all of the cells align well if we split the plotting between conditions, which we can do by adding the split.by argument to the DimPlot() function:
DimPlot(data_integrated, split.by = "orig.ident", group.by = "orig.ident")
ggsave(paste("3_Integration_Figures/Integration_UMAP_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)

####################################################################
## 4. Clustering
####################################################################
# To overcome the extensive technical noise in the expression of any single gene for scRNA-seq data, Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes, with each PC essentially representing a “metagene”? that combines information across a correlated gene set. Determining how many PCs to include in the clustering step is therefore important to ensure that we are capturing the majority of the variation, or cell types, present in our dataset.

# The elbow plot is a helpful way to determine how many PCs to use for clustering so that we are capturing majority of the variation in the data. The elbow plot visualizes the standard deviation of each PC, and we are looking for where the standard deviations begins to plateau. Essentially, where the elbow appears is usually the threshold for identifying the majority of the variation. However, this method can be quite subjective.

# Let’s draw an elbow plot using the top 60 PCs:
# Plot the elbow plot
ElbowPlot(object = data_integrated, ndims = 60)
ggsave(paste("3_Integration_Figures/ElbowPlot_", project_name, ".png", sep = ""), width = 10 , height = 8)

##################################################
## Cluster the cells
##################################################
# Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. A nice in-depth description of clustering methods is provided in the SVI Bioinformatics and Cellular Genomics Lab course.

# We will use the FindClusters() function to perform the graph-based clustering. The resolution is an important argument that sets the “granularity” of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.

# Determine the K-nearest neighbor graph
data_integrated <- FindNeighbors(object = data_integrated, dims = 1:max_dim)

# Determine the clusters for various resolutions                                
data_integrated <- FindClusters(object = data_integrated, resolution = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4))
# If we look at the metadata of our Seurat object(seurat_integrated@meta.data), there is a separate column for each of the different resolutions calculated.

# Explore resolutions only works in Rstudio local not in Nextera or Indar
#data_integrated@meta.data %>% 
#  View()

snn_dim <- rep(NA,11)
snn_dim[1] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.4))[1]
snn_dim[2] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.5))[1]
snn_dim[3] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.6))[1]
snn_dim[4] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.7))[1]
snn_dim[5] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.8))[1]
snn_dim[6] <- dim(table(data_integrated@meta.data$integrated_snn_res.0.9))[1]
snn_dim[7] <- dim(table(data_integrated@meta.data$integrated_snn_res.1))[1]
snn_dim[8] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.1))[1]
snn_dim[9] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.2))[1]
snn_dim[10] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.3))[1]
snn_dim[11] <- dim(table(data_integrated@meta.data$integrated_snn_res.1.4))[1]
names(snn_dim) <- seq(0.4,1.4,0.1)
print(snn_dim)

# To choose a resolution to start with, we often pick something in the middle of the range like 0.6 or 0.8. We will start with a resolution of 0.8 by assigning the identity of the clusters using the Idents() function.
# Assign identity of clusters
resolution_find_clusters <- "integrated_snn_res.0.6"
Idents(object = data_integrated) <- resolution_find_clusters

# choose how many clusters you want and assign it to seurat_clusters
print("WARNING: Choose the resolution that you want.")
data_integrated$seurat_clusters <- data_integrated$integrated_snn_res.0.6

# Plot the UMAP
DimPlot(data_integrated, reduction = "umap", label = TRUE, label.size = 6)
ggsave(paste("3_Integration_Figures/Clusters_UMAP_", project_name, ".png", sep = ""), width = 10 , height = 8)

##################################################
## Segregation of clusters by sample
##################################################
# Explore the distribution of cells per cluster in each sample:
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(data_integrated, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

print(n_cells)
# View table
#View(n_cells)

#Let’s plot the distribution among clusters using our custom function:
png(paste("3_Integration_Figures/Distribution_Cell_", project_name, ".png", sep = ""), width = 700 , height = 700)
plot_integrated_clusters(data_integrated) 
dev.off()

# We can visualize the cells per cluster for each sample using the UMAP: UMAP of cells in each cluster by sample
DimPlot(data_integrated, label = TRUE, split.by = "orig.ident")  + NoLegend()
ggsave(paste("3_Integration_Figures/Clusters_UMAP_bySample_", project_name, ".png", sep = ""), width = 30 , height = 20)
# Generally, we expect to see the majority of the cell type clusters to be present in all conditions; however, depending on the experiment we might expect to see some condition-specific cell types present. These clusters look pretty similar between conditions, which is good since we expected similar cell types to be present in both control and stimulated conditions.

##################################################
## Segregation of clusters by various sources of uninteresting variation
##################################################
# Explore additional metrics, such as the number of UMIs and genes per cell and mitochondrial gene expression by UMAP. 
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rb")

FeaturePlot(data_integrated, reduction = "umap", features = metrics, pt.size = 0.4, order = TRUE, min.cutoff = 'q10', label = TRUE)
ggsave(paste("3_Integration_Figures/Metadata_FeaturePlot_", project_name, ".png", sep = ""), width = 12 , height = 8)

##################################################
## Exploration of the PCs driving the different clusters
##################################################
# Explore how well the clusters separate by the different PCs. To visualize this information, we need to extract the UMAP coordinate information for the cells along with their corresponding scores for each of the PCs to view by UMAP.
# First, we identify the information we would like to extract from the Seurat object, then, we can use the FetchData() function to extract it.

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:max_dim), "ident", "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(data_integrated, vars = columns)

# In the UMAP plots, the cells are colored by their PC score for each respective principal component.
# Let’s take a quick look at the top PCs:
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(data_integrated, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
# NOTE: choose how many PCs you want to analyse 
map(paste0("PC_", 1:max_dim), function(pc){
  ggplot(pc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), alpha = 0.7) +
    scale_color_gradient(guide = FALSE, low = "grey90",  high = "blue")  +
    geom_text(data=umap_label, aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
ggsave(paste("3_Integration_Figures/PC_cluster_relation_", project_name, ".png", sep = ""), width = 12 , height = 10)

# Save integrated seurat object
saveRDS(data_integrated, file = paste("1_SeuratObjects/", "Data_integrated_", project_name, ".rds", sep = ""))

# Write a Log output
Date_time <- Sys.time()
log_data <- as.data.frame(Date_time, row.names = NULL)
log_data$Project_name <- project_name
log_data$Samples <- do.call(paste, c(as.list(samples_v), sep = ", "))
log_data$Organism <- organism
log_data$Approach <- approach
log_data$Integration_reduction <- integration_reduction
log_data$Filtering_parameters <- filterting
log_data$PC_number_doublet <- pc_number_doublet
log_data$PC_number_integration <- max_dim
log_data$Resolution_find_clusters <- resolution_find_clusters

write.table(log_data, file = "2_scPipeline_Integration_13052022_LOG.log", sep = "\t", row.names = FALSE)
