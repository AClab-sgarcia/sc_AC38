## Date: 13/05/2022

## In this pipeline this steps are done: 
##  1. Automaic annotation with SingleR 
##  2. Automaic annotation with scCATCH
##  3. Manual annotation
##  4. Signature annotation

## Directory structure:
##  1_SeuratObject 
##  2_QC_Figures
##  3_Integration_Figures
##  4_MarkerIdentification_Figures

dir.create("4_MarkerIdentification_Figures")
dir.create("4_MarkerIdentification_Figures/SingleR")
dir.create("4_MarkerIdentification_Figures/scCatch")
dir.create("4_MarkerIdentification_Figures/Manual")
dir.create("4_MarkerIdentification_Figures/Signatures")

## General pipeline for single-cell data analysis
## Following:
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Guided%20Clustering%20Tutorial%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Introduction%20to%20scRNA-seq%20integration%20(May2022).R 
## -  Single-cell RNA-seq data analysis workshop: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html 

## Definitions of the functions from: https://cloud.r-project.org/web/packages/Seurat/Seurat.pdf (Version 4.1.1 || Date 2022-05-01)

####################################################################
## Libraries
####################################################################
#  BiocManager::install("SingleR")
suppressPackageStartupMessages(require("SingleR"))
# install.packages("XLConnect")
suppressPackageStartupMessages(require("XLConnect"))
# install.packages("scCATCH")
suppressPackageStartupMessages(require("scCATCH"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("future"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("dplyr"))


# suppressPackageStartupMessages(require("dplyr"))
# suppressPackageStartupMessages(require("scater"))
# suppressPackageStartupMessages(require("patchwork"))
# suppressPackageStartupMessages(require("celldex"))
# suppressPackageStartupMessages(require("purrr"))
# suppressPackageStartupMessages(require("scales"))
# suppressPackageStartupMessages(require("biomaRt"))
# suppressPackageStartupMessages(require("SeuratWrappers"))
# suppressPackageStartupMessages(require("biomaRt"))
# suppressPackageStartupMessages(require("SeuratDisk"))
# suppressPackageStartupMessages(require("celldex"))
# # BiocManager::install("celldex")
# # remotes::install_github("ZJUFanLab/scCATCH")


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

data_integrated <- readRDS(paste("1_SeuratObjects/Data_integrated_", project_name,  ".rds", sep = ""))

if (approach == "Norm_Feature_Scale"){
  DefaultAssay(data_integrated) <- "RNA"
} else if (approach == "SCT") {
  DefaultAssay(data_integrated) <- "SCT"
  data_integrated <- PrepSCTFindMarkers(data_integrated)
  saveRDS(data_integrated, paste("1_SeuratObjects/Data_integrated_PrepSCT_", project_name,  ".rds", sep = ""))
}

###################################################################
## 1. Automatic annotation with SingleR
###################################################################
print("For annotating celldex package its used; HumanPrimaryCellAtlasData for human and MouseRNAseqData for mouse. Check http://bioconductor.org/packages/release/data/experiment/manuals/celldex/man/celldex.pdf for more references.")

# Let’s get reference datasets from celldex package. Note that there are two cell type assignments, label.main and label.fine. We’re only going to run the annotation against the MouseRNAseqData. 
# if (organism == "human") {
#   ref <- celldex::HumanPrimaryCellAtlasData()
# } else if (organism == "mouse"){
#   ref <- celldex::MouseRNAseqData()
# }

if (organism == "human") {
  ref <- readRDS("utils/ref_human.rds")
} else if (organism == "mouse"){
  ref <- readRDS("utils/ref_mouse.rds")
}

data_integrated_automatic <- data_integrated

# Let’s convert our Seurat object to single cell experiment (SCE) for convenience. After this, using SingleR becomes very easy:
sce <- as.SingleCellExperiment(DietSeurat(data_integrated_automatic))

# The finer cell types annotations are you after, the harder they are to get reliably. This is where comparing many databases, as well as using individual markers from literature, would all be very valuable.
main <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.main)
fine <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.fine)

# Lets see the summary of general cell type annotations.
table(main$pruned.labels)
table(fine$pruned.labels)

# Lets add the annotations to the Seurat object metadata so we can use them:
data_integrated_automatic@meta.data$main <- main$pruned.labels
data_integrated_automatic@meta.data$fine <- fine$pruned.labels

# Finally, lets visualize the fine-grained annotations.
data_integrated_automatic_main <- SetIdent(data_integrated_automatic, value = "main")
data_integrated_automatic_fine <- SetIdent(data_integrated_automatic, value = "fine")

DimPlot(data_integrated_automatic_main, label = T , repel = T, label.size = 5) + NoLegend()
ggsave(file= paste("4_MarkerIdentification_Figures/SingleR/AutomaticMarker_SingleR_Main_", project_name, ".png", sep = ""), width = 8, height = 6)

DimPlot(data_integrated_automatic_fine, label = T , repel = T, label.size = 5) + NoLegend()
ggsave(file= paste("4_MarkerIdentification_Figures/SingleR/AutomaticMarker_SingleR_Fine_", project_name, ".png", sep = ""), width = 8, height = 6)

saveRDS(data_integrated_automatic_main, paste("1_SeuratObjects/Annotation_SingleR_Main_", project_name,  ".rds", sep = ""))
saveRDS(data_integrated_automatic_fine, paste("1_SeuratObjects/Annotation_SingleR_Fine_", project_name,  ".rds", sep = ""))

##################################################
## SingleR results QC - https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#1_Introduction
##################################################
# plotScoreHeatmap() displays the scores for all cells across all reference labels, which allows users to inspect the confidence of the predicted labels across the dataset. Ideally, each cell (i.e., column of the heatmap) should have one score that is obviously larger than the rest, indicating that it is unambiguously assigned to a single label. A spread of similar scores for a given cell indicates that the assignment is uncertain, though this may be acceptable if the uncertainty is distributed across similar cell types that cannot be easily resolved.
png(filename = paste("4_MarkerIdentification_Figures/SingleR/QC_AutomaticMarker_SingleR_Main_Heatmap_", project_name, ".png", sep = ""), width = 550, height = 550)
plotScoreHeatmap(main)
dev.off()

png(filename = paste("4_MarkerIdentification_Figures/SingleR/QC_AutomaticMarker_SingleR_Fine_Heatmap_", project_name, ".png", sep = ""), width = 550, height = 550)
plotScoreHeatmap(fine)
dev.off()

# Another diagnostic is based on the per-cell “deltas”, i.e., the difference between the score for the assigned label and the median across all labels for each cell. Low deltas indicate that the assignment is uncertain, which is especially relevant if the cell’s true label does not exist in the reference. 
png(filename = paste("4_MarkerIdentification_Figures/SingleR/QC_AutomaticMarker_SingleR_Main_Distribution_", project_name, ".png", sep = ""), width = 750, height = 750)
plotDeltaDistribution(main, ncol = 3)
dev.off()

png(filename = paste("4_MarkerIdentification_Figures/SingleR/QC_AutomaticMarker_SingleR_Fine_Distribution_", project_name, ".png", sep = ""), width = 750, height = 750)
plotDeltaDistribution(fine, ncol = 3)
dev.off()

# The pruneScores() function will remove potentially poor-quality or ambiguous assignments based on the deltas. The minimum threshold on the deltas is defined using an outlier-based approach that accounts for differences in the scale of the correlations in various contexts - see ?pruneScores for more details. SingleR() will also report the pruned scores automatically in the pruned.labels field where low-quality assignments are replaced with NA.
summary(is.na(main$pruned.labels))
summary(is.na(fine$pruned.labels))

###################################################################
## 2. Automatic annotation with scCATCH - https://cran.r-project.org/web/packages/scCATCH/vignettes/tutorial.html
###################################################################
# For scRNA-seq data, we suggest to revise the gene symbols with rev_gene(). geneinfo is the system data.frame containing the information of human and mouse from NCBI gene (updated in Jan. 2, 2022). To use your own geneinfo data.frame, please refer to demo_geneinfo to build a new one, e.g., rat, zebrafish, Drosophila, C. elegans, etc.
if (organism == "human"){
  species <- "Human"
} else if (organism == "mouse"){
  species <- "Mouse"
}

data_integrated_automatic <- data_integrated

# Revise gene symbols
ref <- rev_gene(data = GetAssayData(data_integrated_automatic, assay = "RNA", slot = "data"), data_type = "data", species = species, geneinfo = geneinfo) # Update Gene symbols in CellMatch according to NCBI Gene symbols (updated in Jan. 2, 2022, https://www.ncbi.nlm.nih.gov/gene).

# Create scCATCH object. Users need to provide the normalized data and the cluster for each cell."
scCATCH_obj <- createscCATCH(data = ref, cluster = as.character(data_integrated_automatic$seurat_clusters))

# Find highly expressed genes with findmarkergene(). Users need to provided the species, tissue, or cancer information. cellmatch is the system data.frame containing the known markers of human and mouse.
write.table(table(cellmatch$tissue), file = "4_MarkerIdentification_Figures/scCatch/scCATCH_TissueTypes.txt", col.names = TRUE)
print("Choose the tissues you are interested in. Tissue list in the following path: 4_MarkerIdentification_Figures/scCatch/scCATCH_TissueTypes.txt")

cellmatch_new <- rbind(cellmatch[cellmatch$species == species & cellmatch$tissue %in% c("Prostate", "Breast", "Mammary epithelium","Blood", "Peripheral blood"), ])
scCATCH_obj <- findmarkergene(object = scCATCH_obj, marker = cellmatch_new, if_use_custom_marker = TRUE)

# Evidence-based score and annotation for each cluster with findcelltype()
scCATCH_obj <- findcelltype(object = scCATCH_obj)
cluster_annotation <- scCATCH_obj@celltype
cluster_annotation <- cluster_annotation[order(cluster_annotation$cluster), ]
write.table(cluster_annotation, file = paste("4_MarkerIdentification_Figures/scCatch/scCATCH_Annotation_", project_name, ".txt", sep = ""), col.names = TRUE, row.names = FALSE)

p <- cluster_annotation %>% 
  ggplot(aes(x=cell_type, fill=cell_type)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Count cell types")
p + coord_flip()
ggsave(paste("4_MarkerIdentification_Figures/scCatch/CellTypes_Barplot_", project_name, ".png", sep = ""), width = 15, height = 8)

for (i in 1:length(data_integrated_automatic@meta.data$seurat_clusters)){
  for (j in 1:length(cluster_annotation$cluster)){
    if (data_integrated_automatic@meta.data$seurat_clusters[i] == cluster_annotation$cluster[j]){
      data_integrated_automatic@meta.data$scScatch_annnotation[i] <- cluster_annotation$cell_type[j]
    }
  }
}

# Finally, lets visualize the fine-grained annotations.
data_integrated_automatic_scCatch <- SetIdent(data_integrated_automatic, value = "scScatch_annnotation")

DimPlot(data_integrated_automatic_scCatch, label = T , repel = T, label.size = 5) + NoLegend()
ggsave(file= paste("4_MarkerIdentification_Figures/scCatch/AutomaticMarker_scCatch_", project_name, ".png", sep = ""), width = 8, height = 6)

saveRDS(data_integrated_automatic_scCatch, paste("1_SeuratObjects/Annotation_scCatch_", project_name,  ".rds", sep = ""))

# remove variables from environment
ls()
rm(data_integrated_automatic)
rm(data_integrated_automatic_fine)
rm(data_integrated_automatic_main)
rm(data_integrated_automatic_scCatch)
rm(scCATCH_obj)
rm(ref)
rm(sce)
rm(main)
rm(fine)
ls()

###################################################################
## 3. Mannual annotation
###################################################################
# For certain functions, each worker needs access to certain global variables. If these are larger than the default limit, you will see this error. To get around this, you can set options(future.globals.maxSize = X), where X is the maximum allowed size in bytes. So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). Note that this will increase your RAM usage so set this number mindfully.
options(future.globals.maxSize = 2000 * 1024^2)
plan("multicore", workers = 18)

data_markers <- FindAllMarkers(data_integrated, min.pct = 0.1, logfc.threshold = 0.25)
dim(data_markers)
table(data_markers$cluster)
top15_markers <- as.data.frame(data_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC))
write.table(top15_markers, file = paste("4_MarkerIdentification_Figures/Manual/Annotation_Manual_top10genesFindAllMarkers_", project_name,  ".txt", sep = ""), row.names = FALSE)

if (organism == "human"){
  features <- c("FLT3", "ITGAM", "CSF1R", "CD14", "CD3E", "CD3D", "CD3G", "CD4", "FOXP3", "IL2RA", "CD8A", "NCR1", "CD79A", "CD19", "S100A8", "S100A9", "CD68", "CD80", "APOE", "EPCAM", "PECAM1", "PDGFRA", "PDGFRB", "THY1", "RGS5", "KRT5", "KRT8" )
} else if(organism == "mouse"){
  features <- c("Flt3", "Itgam", "Csf1r", "Cd14", "Cd3e", "Cd3d", "Cd3g", "Cd4", "Foxp3", "Il2ra", "Cd8a", "Ncr1", "Cd79a", "Cd19", "S100a8", "S100a9", "Cd68", "Cd80", "Apoe", "Epcam", "Pecam1", "Pdgfra", "Pdgfrb", "Thy1", "Rgs5", "Krt5", "Krt8")
}

# In the features line write the markers for the Cell Type
for (feature in features){
  # png(paste("4_MarkerIdentification_Figures/Manual/Annotation_Manual_FeaturePlot_", feature , "_", project_name, ".png", sep = ""), width = 8, height = 6)
  FeaturePlot(data_integrated, 
              reduction = "umap", 
              features = feature, 
              order = TRUE,
              min.cutoff = 'q10', 
              label = TRUE,
              keep.scale = "all")
  # dev.off()
  
  ggsave(file= paste("4_MarkerIdentification_Figures/Manual/Annotation_Manual_FeaturePlot_", feature , "_", project_name, ".png", sep = ""), width = 8, height = 6)
}

print("WARNIING: Make the changes in the next variable for the cluster names.")

new.cluster.ids <- c("0" = "0_",
                     "1" = "1_", 
                     "2" = "2_",
                     "3" = "3_",
                     "4" = "4_",
                     "5" = "5_",
                     "6" = "6_",
                     "7" = "7_",
                     "8" = "8_",
                     "9" = "9_",
                     "10" = "10_",
                     "11" = "11_",
                     "12" = "12_",
                     "13" = "13_",
                     "14" = "14_",
                     "15" = "15_",
                     "16" = "16_",
                     "17" = "17_", 
                     "18" = "18_", 
                     "19" = "19_", 
                     "20" = "20_",
                     "21" = "21_",
                     "22" = "22_",
                     "23" = "23_",
                     "24" = "24_",
                     "25" = "25_")

names(new.cluster.ids) <- levels(data_integrated)
data_annotated <- RenameIdents(data_integrated, new.cluster.ids)

png("4_MarkerIdentification_Figures/Manual/Manual annotation.png", width = 1500, height = 1000)
DimPlot(data_annotated, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.8) 
dev.off()

##################################################################################################
# my own classification of clusters
##################################################################################################
if (organism == "human"){
  T_Cells <- c("CD3E", "CD3D", "CD3G", "CD4", "CD8A") 
  NK <- c("NCR1")
  B_cell <- c("CD79A", "CD19")
  Neutrophils <- c("S100A8", "S100A9")
  Macrophage <- c("CD68")
  Epithelial <- c("EPCAM")
  Endothelial <- c("PECAM1", "KRT8", "KRT5")
  Fibroblast <- c("THY1", "PDGFRB", "PDGFRA", "THY1") 
  Mural_cells <- c("RGS5")
    Dendritic <- c("FLT3")
} else if(organism == "mouse"){
  T_Cells <- c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a") 
  NK <- c("Ncr1")
  B_cell <- c("Cd79a", "Cd19")
  Neutrophils <- c("S100a8", "S100a9")
  Macrophage <- c("Cd68", "Igtam", "Csf1r")
  Epithelial <- c("Epcam")
  Endothelial <- c("Pecam1", "Krt8", "Krt5")
  Fibroblast <- c("Thy1", "Pdgfrb", "Pdgfra", "Thy1") 
  Mural_cells <- c("Rgs5")
    Dendritic <- c("Flt3")
}

# plot ALL markers in dotplot
to_plot <- unique(c(T_Cells, Epithelial, Fibroblast, NK, Mural_cells, Endothelial, Macrophage, B_cell, Neutrophils))
to_plot <- to_plot[which(to_plot %in% row.names(data_annotated@assays$SCT))]
data_annotated@meta.data$seurat_clusters <- data_annotated$seurat_clusters

pdf(paste("4_MarkerIdentification_Figures/Manual/DotPlot_Many_markers_", project_name, ".pdf",sep=''),width=15,height=15)
DotPlot(data_annotated, features = to_plot) 
dev.off()

# Violin plots can also be split on some variable. Simply add the splitting variable to object metadata and pass it to the split.by argument
pdf(paste("4_MarkerIdentification_Figures/Manual/VlnPlot_Many_markers_", project_name, ".pdf",sep=''),width=28,height=28)
VlnPlot(data_annotated, features = to_plot)
dev.off()

# Single cell heatmap of feature expression
pdf(paste("4_MarkerIdentification_Figures/Manual/Heatmap_Many_markers_", project_name, ".pdf",sep=''),width=12,height=12)
DoHeatmap(subset(data_annotated, downsample = 100), features = to_plot, size = 3)
dev.off()

saveRDS(data_annotated, paste("1_SeuratObjects/Annotation_Manual_", project_name,  ".rds", sep = ""))

## LOCAL
###################################################################
## 4. Signature annotation
###################################################################
data_integrated <- readRDS(paste("x:/sgarcia/sc_AC38_AC50/1_SeuratObjects/Data_integrated_PrepSCT_", project_name,  ".rds", sep = ""))

if (organism == "human"){
  wb <- loadWorkbook("Bibliography/Human/Human_Markers.xlsx") 
  markers <- readWorksheet(wb, sheet = getSheets(wb))
} else if (organism == "mouse") {
  wb <- loadWorkbook("Bibliography/Mouse/Mouse_Markers.xlsx") 
  markers <- readWorksheet(wb, sheet = getSheets(wb))
  saveWorkbook(wb)
}

print("Check the README dataframe in 'markers' variable to choose which signatures to use. Write the name of the sheet in the 'signatures' variable.")
print(markers$README)
sheets <- c("1", "2", "3", "4", "5", "6")

for (sheet in sheets){
  
  print(sheet)
  marker_df <- markers[[sheet]]
  
  graphs_signature <- list()
  genes <- rownames(data_integrated@assays$SCT@counts)
  
  for(i in 1:ncol(marker_df)){
    
    signature <- marker_df[ ,i, drop=FALSE]
    
    function_name <- colnames(signature)
    graph_name <- paste(function_name,"1",sep="")
    
    signature <- signature[[function_name]]
    # signature <- tools::toTitleCase(tolower(signature))
    signature <- intersect(genes,signature)
    
    data <- AddModuleScore(object = data_integrated,
                           features = list(signature),
                           name = function_name)
    
    X <- VlnPlot(data, features = graph_name, sort = T)
    Y <- FeaturePlot(data, features=graph_name, min.cutoff = "q10", pt.size = 0.3)
    Z <- plot_grid(X, Y, labels = "AUTO")
    
    graphs_signature[[i]] <- Z
    names(graphs_signature)[i] <- function_name
    
  }
  
  pdf(paste("X:/sgarcia/sc_AC38_AC50/4_MarkerIdentification_Figures/Signatures/", sheet, ".pdf", sep = ""), height = 6, width = 10)
  print(graphs_signature)
  dev.off()

}





