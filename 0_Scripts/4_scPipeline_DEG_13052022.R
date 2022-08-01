## Date: 13/05/2022

## In this pipeline this steps are done: 
##  1. Differential expression testing - among conditions
##  2. Differential expression testing - among conditions MAST
##  3. Differential expression testing - pseudobulk
##  4. CNV inference

## Directory structure:
##  5_DEG

dir.create("5_DEG")
dir.create("5_DEG/Wilcox")
dir.create("5_DEG/MAST")
dir.create("5_DEG/Pseudobulk")

## General pipeline for single-cell data analysis
## Following:
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Guided%20Clustering%20Tutorial%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Introduction%20to%20scRNA-seq%20integration%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Differential%20expression%20testing%20(May2022).R
## -  Single-cell RNA-seq data analysis workshop: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html 

## Definitions of the functions from: https://cloud.r-project.org/web/packages/Seurat/Seurat.pdf (Version 4.1.1 || Date 2022-05-01)

####################################################################
## Libraries
####################################################################
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("MAST"))
suppressPackageStartupMessages(require("DESeq2"))
suppressPackageStartupMessages(require("scater"))

# suppressPackageStartupMessages(require("dplyr"))
# suppressPackageStartupMessages(require("patchwork"))
# suppressPackageStartupMessages(require("cowplot"))
# suppressPackageStartupMessages(require("ggplot2"))
# suppressPackageStartupMessages(require("celldex"))
# suppressPackageStartupMessages(require("purrr"))
# suppressPackageStartupMessages(require("SingleR"))
# suppressPackageStartupMessages(require("scales"))
# suppressPackageStartupMessages(require("biomaRt"))

###################################################################
## Variable prep
###################################################################
# The variables are stored in the log file from the 1_scPipeline_QC
log_file <- read.delim("2_scPipeline_Integration_13052022_LOG.log", header = TRUE)

project_name <- log_file$Project_name
samples_v <- log_file$Samples
organism <- log_file$Organism
approach <- log_file$Approach
integration_reduction <- log_file$Integration_reduction
filterting <- log_file$Filtering_parameters
pc_number_doublet <- log_file$PC_number_doublet
pc_number_integration <- log_file$PC_number_integration
resolution <- log_file$Resolution_find_clusters

print("Several annotation types have been done in the previous step, choose which one to use between: scCatch, SingleR_Main, SingleR_Fine, Manual")

data_annotated <- readRDS(paste("1_SeuratObjects/Annotation_Manual_", project_name,  ".rds", sep = ""))

# Add a pheno column to metadata with the info of the contrasts for DEGs. 
data_annotated[["pheno"]] <- NA
data_annotated$pheno[which(data_annotated$orig.ident == "S_01")] <- "PtenHo"
data_annotated$pheno[which(data_annotated$orig.ident == "S_02")] <- "PtenHo_Ppargc1aHo"
data_annotated$pheno[grep("HoHo", data_annotated$orig.ident)] <- "PtenHo_Stk11Ho"
data_annotated$pheno[grep("HoHe", data_annotated$orig.ident)] <- "PtenHo_Stk11He"
data_annotated$pheno[grep("lung$", data_annotated$orig.ident)] <- "PtenHo_Stk11Ho_lung"

table(data_annotated$orig.ident)
table(data_annotated$pheno)

###################################################################
## 1. Differential expression testing - among conditions Wilcox
###################################################################
min_diff_scdegs <- 0; # min % of cells per cluster expressing the gene in
total_clusters <- as.character(sort(unique(data_annotated$seurat_clusters)))
combinations <- combn(unique(data_annotated$pheno), 2)

for (comb in 1:dim(combinations)[2]){ 
  sc_deg_df <- c()
    for(clus in 1:length(total_clusters)){
      cluster_subset <- subset(data_annotated, subset=seurat_clusters == total_clusters[clus])
      Idents(cluster_subset) <-cluster_subset$pheno
      #DEGS both up and down
      if (approach == "SCT"){
        tryCatch( 
          cluster_subset <- PrepSCTFindMarkers(cluster_subset), error = function(e) return(NA)
          )
      } else {
        cluster_subset <- cluster_subset
      }
      tryCatch(
        met_markers <- FindMarkers(cluster_subset, ident.1 = combinations[1,comb], ident.2 = combinations[2,comb], verbose = FALSE, only.pos= FALSE), error = function(e) return(NA)
      )
      
      met_markers <- met_markers[(met_markers$p_val_adj <= 0.05) & (abs(met_markers$pct.2-met_markers$pct.1) >= min_diff_scdegs),]
      sc_deg_df <- rbind(sc_deg_df, cbind(rep(total_clusters[clus], dim(met_markers)[1]), row.names(met_markers), met_markers))
    }
  colnames(sc_deg_df)[1:2] <- c("Seurat_clusters","Genes")
  write.table(sc_deg_df, paste("5_DEG/Wilcox/Wilcox_", combinations[1, comb],"_vs_", combinations[2, comb], ".txt", sep=''), row.names = FALSE)
  assign(paste(combinations[1,comb],"_",combinations[2,comb],"_DEG_df", sep=''), sc_deg_df)
}

###################################################################
## 2. Differential expression testing - among conditions MAST
###################################################################
min_diff_scdegs <- 0; # min % of cells per cluster expressing the gene in
total_clusters <- as.character(sort(unique(data_annotated$seurat_clusters)))
combinations <- combn(unique(data_annotated$pheno), 2)

for (comb in 1:dim(combinations)[2]){ 
  sc_deg_df <- c()
  for(clus in 1:length(total_clusters)){
    cluster_subset <- subset(data_annotated, subset=seurat_clusters == total_clusters[clus])
    Idents(cluster_subset) <-cluster_subset$pheno
    #DEGS both up and down
    if (approach == "SCT"){
      tryCatch( 
        cluster_subset <- PrepSCTFindMarkers(cluster_subset), error = function(e) return(NA)
      )
    } else {
      cluster_subset <- cluster_subset
    }
    tryCatch(
      met_markers <- FindMarkers(cluster_subset, ident.1 = combinations[1,comb], ident.2 = combinations[2,comb], test.use = "MAST", verbose = FALSE, only.pos= FALSE), error = function(e) return(NA)
    )
    
    met_markers <- met_markers[(met_markers$p_val_adj <= 0.05) & (abs(met_markers$pct.2-met_markers$pct.1) >= min_diff_scdegs),]
    sc_deg_df <- rbind(sc_deg_df, cbind(rep(total_clusters[clus], dim(met_markers)[1]), row.names(met_markers), met_markers))
  }
  colnames(sc_deg_df)[1:2] <- c("Seurat_clusters","Genes")
  write.table(sc_deg_df, paste("5_DEG/MAST/MAST_", combinations[1, comb],"_vs_", combinations[2, comb], ".txt", sep=''), row.names = FALSE)
  assign(paste(combinations[1,comb],"_",combinations[2,comb],"_DEG_df", sep=''), sc_deg_df)
}

###################################################################
## 3. Differential expression testing - pseudobulk
###################################################################
sce <- as.SingleCellExperiment(DietSeurat(data_annotated))

# The scater package provides a variey of tools for preprocessing and quality control of single-cell transcriptomic data. For completeness, we will apply some minimal filtering steps to
#   - remove undetected genes
#   - remove cells with very few or many detected genes
#   - remove very lowly expressed genes
#   - compute normalized expression values for visualization

# For more thorough preprocessing, we refer to the Quality control with scater vignette.
# remove undetected genes

# no es necesario ya que ya tenemos esa info 
dim(sce)
sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# We use calculateQCMetrics (it's defunct, use addPerCellQC instead) to compute various quality control metrics for each cell and gene, stored in the colData and rowData, respectively, and proceed with filtering cells and genes as noted above:
# calculate quality control (QC) metrics
sce <- addPerCellQC(sce)


###################################################################
## 4. CNV inference
###################################################################

#### CNV inference http://www.bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix="../inst/extdata/oligodendroglioma_expression_downsampled.counts.matrix.gz",
  annotations_file="../inst/extdata/oligodendroglioma_annotations_downsampled.txt",
  delim="\t",
  gene_order_file="../inst/extdata/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
  ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

out_dir = tempfile()
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out_dir,
  cluster_by_groups=TRUE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  png_res=60
)

