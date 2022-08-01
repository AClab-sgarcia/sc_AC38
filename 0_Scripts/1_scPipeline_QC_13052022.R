## Date: 13/05/2022

## In this pipeline this steps are done: 
##  1. Create Seurat object
##  2. Merging samples
##  3. QC analysis (plots)
##  4. Filtering
##  5. Gene filtering
##  6. Cell-cycle scoring
##  7. Doublet Finder and the automatic annotation of each sample

## Directory structure:
##  1_SeuratObject 
##  2_QC_Figures

dir.create("1_SeuratObjects")
dir.create("2_QC_Figures")
dir.create("2_QC_Figures/DoubletFinder")

## General pipeline for single-cell data analysis
## Following:
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Guided%20Clustering%20Tutorial%20(May2022).R 
## -  https://github.com/AClab-sgarcia/Single-cell/blob/main/Seurat%20-%20Introduction%20to%20scRNA-seq%20integration%20(May2022).R 
## -  Single-cell RNA-seq data analysis workshop: https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html 
## -  https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html 

## Definitions of the functions from: https://cloud.r-project.org/web/packages/Seurat/Seurat.pdf (Version 4.1.1 || Date 2022-05-01)

####################################################################
## Libraries
####################################################################
suppressPackageStartupMessages(require("dplyr"))
# install.packages("scater")
suppressPackageStartupMessages(require("scater"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("ggplot2"))
# BiocManager::install("celldex")
suppressPackageStartupMessages(require("celldex"))
suppressPackageStartupMessages(require("purrr"))
# BiocManager::install("SingleR")
suppressPackageStartupMessages(require("SingleR"))
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
suppressPackageStartupMessages(require("DoubletFinder"))
suppressPackageStartupMessages(require("scales"))
# BiocManager::install("biomaRt")
# suppressPackageStartupMessages(require("biomaRt"))
suppressPackageStartupMessages(require("data.table"))

###################################################################
## Variable prep
###################################################################
# Big data
input_path <- file.path("X:/sgarcia/sc_AC38_AC50/RAW data")

# for Nextera (inside the mounted docker)
input_path <- file.path("RAW data")

# Cell ranger filtered matrix access 
path_matrix <- ("outs/filtered_feature_bc_matrix/")

# Organism (options: "human" || "mouse")
organism <- "mouse"

# Project name
project_name <- "AC38"

# for AC38
path_cellRanger <- paste(input_path, "AC-38_10x-3RNA/CellRangerCount", sep="/")
# Sample names
samples_v <- c("S_01", "S_02")

# Approach (options: "Norm_Feature_Scale" || "SCT")
approach <- "SCT"
print("In this script SCT only regresses out mithocondrial percentage, change accordingly to your needs :)")

# Filtering paramenters (use after QC)
# min_nFeature == nGene
# UMI == unique RNA molecules
min_nFeature <- 200
min_log10GenesPerUMI <- 0.8
min_percent_mt <- 20

####################################################################
## 1. Analyzing CellRangerCount's web_summar for the Sample
####################################################################
## Summary:
##  - Include intons: FALSE
##  - Reference transcriptome: /vols/GPArkaitz_bigdata/DATA_shared/Genomes/refdata-gex-mm10-2020-A
##  - Pipeline Version: cellranger-6.1.0

## Number of cells per sample according to Cell Ranger
## S_01   S_02   
## 8.573  8.553

####################################################################
## 2. Set-up the Seurat Object
####################################################################
for (sample in samples_v){
  path_to_sample <- paste(path_cellRanger, sample, path_matrix, sep = "/")
  
  # Read10X: Enables easy loading of sparse data matrices provided by 10X genomics. Returns a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
  # the data.dir variable: Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X.
  seurat_data <- Read10X(data.dir = path_to_sample) 

  # Initialize the Seurat object with the raw (non-normalized data). The object serves as a container that contains both data (like the count matrix) and  analysis (like PCA, or clustering results) for a single-cell dataset.
  # - min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well.
  # - min.features: Include cells where at least this many features are detected.
  # (https://www.biostars.org/p/407339/) min.cells helps limit the number of genes used by removing those unlikely to play any part in differentiating groups of cells due to being expressed in very few cells. In general, most genes removed will be those with zero counts across all cells. min.features removes dead cells and empty droplets where few genes are detected.
  seurat_obj <- CreateSeuratObject(counts = seurat_data, project = sample, min.cells = 0, min.features = 0)
  assign(sample, seurat_obj)
}

for (sample in samples_v){
  # We are going to save each of the objects separately 
  saveRDS(get(sample), file = paste("1_SeuratObjects/", sample, "_", project_name, ".rds", sep = ""))
}

## If you are working with raw data matrixes use the following chunk:
# project_seurat_obj <- CreateSeuratObject(counts = data_raw, min.cells = 0, min.features = min_nFeature, project =project_name)

####################################################################
## 3. Merging Seurat Objects - https://satijalab.org/seurat/articles/merge_vignette.html
####################################################################
# merge() merges the raw count matrices of two or mores Seurat objects and creates a new Seurat object with the resulting combined raw count matrix. To easily tell which original object any particular cell came from, you can set the add.cell.ids parameter with an c(x, y) vector.
# for more than 2 objects: y = c(sample2, sample3, ...)
if (length(samples_v) > 1){
  rawData_combined <- merge(x = get(samples_v[1]), y = sapply(samples_v[-1], get), add.cell.ids = samples_v, project = project_name)
}

# Explore merged metadata
# View(rawData_combined@meta.data)

# Notice the cell names now have an added identifier
head(colnames(rawData_combined))

# Number of cells in each sample, should be the same as in the CellRanger unless min.cells or min.features has been used in the CreateSeuratObject function. 
table(rawData_combined$orig.ident)

####################################################################
## 4. Standard pre-processing workflow - QC and filtering
####################################################################
##################################################
## QC and selecting cells for further analysis
##################################################
# Mitochondrial Ratio: calculate and save it in the percent_mt column of the data (percentage of reads that map to the mitochondrial genome). 
# As for human and mouse it's different, let's make it general with grep function over the rownames of the data and use ignore.case = TRUE.
if (organism == "human") {
  rawData_combined[["percent_mt"]] <- PercentageFeatureSet(rawData_combined, pattern = "^MT-")
} else if (organism == "mouse"){
  rawData_combined[["percent_mt"]] <- PercentageFeatureSet(rawData_combined, pattern = "^mt-")
}

# Ribosomal Ratio
if (organism == "human") {
  rawData_combined[["percent_rb"]] <- PercentageFeatureSet(rawData_combined, pattern = "^RP[SL]")
} else if (organism == "mouse"){
  rawData_combined[["percent_rb"]] <- PercentageFeatureSet(rawData_combined, pattern = "^Rp[sl]")
}

# Add number of genes per UMI for each cell to metadata
rawData_combined$log10GenesPerUMI <- log10(rawData_combined$nFeature_RNA) / log10(rawData_combined$nCount_RNA)

VlnPlot(rawData_combined, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_rb"), ncol = 4)
ggsave(paste("2_QC_Figures/QC_VlnPlot_", project_name, ".png", sep = ""), width = 16, height = 10)
# d + geom_jitter(shape = 13, alpha = 0.2, aes(colour = rep("grey", length(rawData_combined@meta.data$orig.ident))))

plot1 <- FeatureScatter(rawData_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(rawData_combined, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot3 <- FeatureScatter(rawData_combined, feature1 = "nCount_RNA", feature2 = "percent_rb")
plot1 + plot2 +plot3
ggsave(paste("2_QC_Figures/QC_FeatureScatter_", project_name, ".png", sep = ""), width = 12 , height = 9)

# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the total number of molecules detected within a cell. Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet). 

# Visualize the number of cell counts per sample (origin identification)
rawData_combined@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  geom_text(stat='count', aes(label=..count..), vjust=1.6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave(paste("2_QC_Figures/QC_NCells_Barplot_", project_name, ".png", sep = ""), width = 10 , height = 8)

# Visualize the number UMIs/transcripts per cell
rawData_combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500, color = "grey") # The UMI counts per cell should generally be above 500, that is the low end of what we expect. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
ggsave(paste("2_QC_Figures/QC_nCountRNA_Density_", project_name, ".png", sep = ""), width = 12 , height = 9)

# Visualize the distribution of genes detected per cell via histogram
rawData_combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300, color = "grey")
ggsave(paste("2_QC_Figures/QC_nFeatureRNA_Density_", project_name, ".png", sep = ""), width = 12 , height = 9)

# Visualize the distribution of genes detected per cell via boxplot
rawData_combined@meta.data %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave(paste("2_QC_Figures/QC_nFeatureRNA_Boxplot_", project_name, ".png", sep = ""), width = 12 , height = 9)

# UMIs vs. genes detected
# Two metrics that are often evaluated together are the number of UMIs and the number of genes detected per cell. We plot the number of genes versus the number of UMIs coloured by the fraction of mitochondrial reads. Mitochondrial read fractions are only high in particularly low count cells with few detected genes (darker colored data points). This could be indicative of damaged/dying cells whose cytoplasmic mRNA has leaked out through a broken membrane, and thus, only mRNA located in the mitochondria is still conserved. These cells are filtered out by our count and gene number thresholds. Jointly visualizing the count and gene thresholds shows the joint filtering effect.
# Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs.
# With this plot we also evaluate the slope of the line, and any scatter of data points in the bottom right hand quadrant of the plot. These cells have a high number of UMIs but only a few number of genes. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
rawData_combined@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent_mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500, color = "grey") +
  geom_hline(yintercept = 250, color = "grey") +
  facet_wrap(~orig.ident)
ggsave(paste("2_QC_Figures/QC_NCountRNA_nFeatureRNA_mt_Scatter_", project_name, ".png", sep = ""), width = 12 , height = 9)

# Mitochondrial counts ratio
# This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark, unless of course you are expecting this in your sample.
# Visualize the distribution of mitochondrial gene expression detected per cell
rawData_combined@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent_mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20, color = "grey")
ggsave(paste("2_QC_Figures/QC_percentmt_Density_", project_name, ".png", sep = ""), width = 12 , height = 9)

# Complexity
# We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts from those lower number of genes over and over again. These low complexity (low novelty) cells could represent a specific cell type (i.e. red blood cells which lack a typical transcriptome), or could be due to some other strange artifact or contamination. Generally, we expect the novelty score to be above 0.80 for good quality cells.

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
rawData_combined@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8, color = "grey")
ggsave(paste("2_QC_Figures/QC_log10GenesPerUMI_Density_", project_name, ".png", sep = ""),  width = 12 , height = 9)

##################################################
## Filtering
##################################################
# Filter according to the plots above
print(paste("WARNING! These filters are applied: nFeature_RNA >", min_nFeature, "& log10GenesPerUMI >", min_log10GenesPerUMI, "& percent_mt <", min_percent_mt, sep = " "))
data_combined <- subset(rawData_combined, subset = nFeature_RNA > min_nFeature & log10GenesPerUMI > min_log10GenesPerUMI & percent_mt < min_percent_mt)

# We can plot the same qc-stats to see the results of the filtering
VlnPlot(data_combined, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_rb"), ncol = 4)
ggsave(paste("2_QC_Figures/QC_VlnPlot_", project_name, "_Filtered.png", sep = ""), width = 16, height = 10)

# We can now see haw many cells we have left
table(Idents(rawData_combined))
table(Idents(data_combined))

# Gene-level filtering
# Within our data we will have many genes with zero counts. These genes can dramatically reduce the average expression for a cell and so we will remove them from our data. We will start by identifying which genes have a zero count in each cell:
# Extract counts
counts <- GetAssayData(object = data_combined, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Now, we will perform some filtering by prevalence. If a gene is only expressed in a handful of cells, it is not particularly meaningful as it still brings down the averages for all other cells it is not expressed in. For our data we choose to keep only genes which are expressed in 10 or more cells. By using this filter, genes which have zero counts in all cells will effectively be removed.
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

data_combined <- data_combined[keep_genes, ]

print(paste("Gene filtering (at least 10 cells) step is removing: ", sum(keep_genes == FALSE), " and ", sum(keep_genes == TRUE), " are mantained.", sep = ""))

# Finally, take those filtered counts and create a new Seurat object for downstream analysis.
# Reassign to filtered Seurat object
split_data <- SplitObject(data_combined, split.by = "orig.ident")

##################################################
## Sample sex
##################################################
# When working with human or animal samples, you should ideally constrain you experiments to a single sex to avoid including sex bias in the conclusions. However this may not always be possible. By looking at reads from chromosomeY (males) and XIST (X-inactive specific transcript) expression (mainly female) it is quite easy to determine per sample which sex it is. It can also bee a good way to detect if there has been any sample mixups, if the sample metadata sex does not agree with the computational predictions.

# To get choromosome information for all genes, you should ideally parse the information from the gtf file that you used in the mapping pipeline as it has the exact same annotation version/gene naming. However, it may not always be available, as in this case where we have downloaded public data. Hence, we will use biomart to fetch chromosome information. As the biomart instances quite often are unresponsive, you can try the code below, but if it fails, we have the file with gene annotations on github here. Make sure you put it at the correct location for the path genes.file to work.

# initialize connection to mart, may take some time if the sites are unresponsive.
# if (organism == "human") {
#   mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# } else if (organism == "mouse"){
#   mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# }

if (organism == "human") {
  mart <- readRDS("utils/mart_human.rds")
} else if (organism == "mouse"){
  mart <- readRDS("utils/mart_mouse.rds")
}

# fetch chromosome info plus some other annotations
# genes_table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart, useCache = F))

if (organism == "human") {
  genes_table <- readRDS("utils/genes_table_human.rds")
} else if (organism == "mouse"){
  genes_table <- readRDS("utils/genes_table_mouse.rds")
}

genes_table <- genes_table[genes_table$external_gene_name %in% rownames(data_combined), ]

# Now that we have the chromosome information, we can calculate per cell the proportion of reads that comes from chromosome Y.
chrY_gene <- genes_table$external_gene_name[genes_table$chromosome_name == "Y"]
data_combined$pct_chrY <- colSums(data_combined@assays$RNA@counts[chrY_gene, ])/colSums(data_combined@assays$RNA@counts)

# Then plot XIST expression vs. chrY proportion. As you can see, the samples are clearly on either side, even if some cells do not have detection of either.
# Xist gene is expressed from the inactivated X chromosomes only, therefore:
# Female: high Xist and low Y
# Male: low Xist and high Y
# Plot as violins.
if (organism == "human") {
  FeatureScatter(data_combined, feature1 = "XIST", feature2 = "pct_chrY")
  ggsave(paste("2_QC_Figures/Sex_FeatureScatter_", project_name, ".png", sep = ""), width = 7 , height = 4)
  VlnPlot(data_combined, features = c("XIST", "pct_chrY"))
  ggsave(paste("2_QC_Figures/Sex_VlnPlot_", project_name, ".png", sep = ""), width = 7 , height = 4)
  
} else if (organism == "mouse"){
  FeatureScatter(data_combined, feature1 = "Xist", feature2 = "pct_chrY")
  ggsave(paste("2_QC_Figures/Sex_FeatureScatter_", project_name, ".png", sep = ""), width = 7 , height = 4)
  VlnPlot(data_combined, features = c("Xist", "pct_chrY"))
  ggsave(paste("2_QC_Figures/Sex_VlnPlot_", project_name, ".png", sep = ""), width = 7 , height = 4)
}

##################################################
## Cell-Cycle Scoring 
##################################################
cc_data <- NormalizeData(data_combined)
cc_data <- CellCycleScoring(object = cc_data, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# View cell cycle scores and phases assigned to cells                                 
# View(cc_data@meta.data)  
# Three columns are added: "S.Score" "G2M.Score" "Phase" 

# Identify the most variable genes
cc_data <- FindVariableFeatures(cc_data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Scale the counts
cc_data <- ScaleData(cc_data)

# Perform PCA
cc_data <- RunPCA(cc_data)

# Plot the PCA colored by cell cycle phase
DimPlot(cc_data, reduction = "pca", group.by= "Phase", split.by = "Phase")
ggsave(paste("2_QC_Figures/CC_DimPlot_", project_name, ".png", sep = ""), width = 12, height = 8)

# We can now plot a violin plot for the cell cycle scores as well.
VlnPlot(cc_data, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", ncol = 2, pt.size = 0.1)
ggsave(paste("2_QC_Figures/CC_VlnPlot_", project_name, ".png", sep = ""), width = 12, height = 8)

# save cell cycle parameters in the meta data
data_combined[[c("S.Score", "G2M.Score", "Phase")]] <- cc_data[[c("S.Score", "G2M.Score", "Phase")]]

saveRDS(data_combined, paste("1_SeuratObjects/QC_data_combined_", project_name,  ".rds", sep = ""))

##################################################
## DoubletFinder
##################################################
# pK Identification (no ground-truth)
# No ground truth means that we dont know if a cell is a single o a doublet
# There are experimental methods to know this, but normally people don't know this
split_data <- SplitObject(data_combined, split.by = "orig.ident")

dim_v <- c()

## for SingleR
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

for (i in 1:length(split_data)) {
  
  if (approach == "Norm_Feature_Scale"){
    x <- NormalizeData(split_data[[i]])
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x)
    sct_param <- FALSE
  } else if (approach == "SCT"){
    x <- SCTransform(split_data[[i]], vars.to.regress = c("percent_mt"))
    sct_param <- TRUE
  }
  
  x <- RunPCA(x)
  ElbowPlot(x, ndims =  60)
  ggsave(paste("2_QC_Figures/DoubletFinder/DF_ElbowPlot_", samples_v[i], "_", project_name, ".png", sep = ""), width = 5 , height = 4)
  
   # From https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
  # Determine percent of variation associated with each PC
  pct <- x[["pca"]]@stdev / sum(x[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim <- min(co1, co2)
  
  dim_v <- c(dim_v, dim)
  x <- FindNeighbors(object = x, dims = 1:dim)
  x <- FindClusters(object = x)
  
  x <- RunUMAP(x, dims = 1:dim, reduction = "pca", metric = "euclidean", min.dist = 0.6)
  sweep_res_x <- paramSweep_v3(x, PCs = 1:dim, sct = sct_param)
  
  sweep_x <- summarizeSweep(sweep_res_x, GT = FALSE)
  bcmvn_x <- find.pK(sweep_x)
  ggplot(bcmvn_x, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  ggsave(paste("2_QC_Figures/DoubletFinder/DF_Doublets_", samples_v[i], "_", project_name, ".png", sep = ""), width = 8 , height = 6)
  
  pK <- bcmvn_x %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- x@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(x@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))
  
  x <- doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = sct_param)
  
  # visualize doublets
  colnames(x@meta.data)
  colname <- colnames(x@meta.data[grep("DF.classifications", colnames(x@meta.data))])
  DimPlot(x, reduction = 'umap', group.by = colname)
  ggsave(file= paste("2_QC_Figures/DoubletFinder/DF_UMAP_", samples_v[i], "_", project_name, ".png", sep = ""), width = 8, height = 6)
  
  # number of singlets and doublets
  table(x@meta.data[colname]) 

  ##################################################
  ## Automaic annotation with SingleR
  ##################################################
  # Let’s convert our Seurat object to single cell experiment (SCE) for convenience. After this, using SingleR becomes very easy:
  sce <- as.SingleCellExperiment(DietSeurat(x))

  main <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref$label.main)

  # Let’s see the summary of general cell type annotations. These match our expectations (and each other) reasonably well.
  table(main$pruned.labels)
  
  # Lets add the annotations to the Seurat object metadata so we can use them:
  x@meta.data$main <- main$pruned.labels

  # Finally, lets visualize the fine-grained annotations.
  x <- SetIdent(x, value = "main")
  
  DimPlot(x, label = T , repel = T, label.size = 5) + NoLegend()
  ggsave(file= paste("2_QC_Figures/DoubletFinder/DF_UMAP_SingleR_", samples_v[i], "_", project_name, ".png", sep = ""), width = 8, height = 6)
  
  saveRDS(x, paste("1_SeuratObjects/QC_", samples_v[i], "_", project_name, ".rds", sep = ""))
  
  ## Filter doublets
  dim(x)
  x <- x[, x@meta.data[, colname] == "Singlet"]
  dim(x)
  saveRDS(x, paste("1_SeuratObjects/QC_Singlet", samples_v[i], "_", project_name,  ".rds", sep = ""))
}

# Write a Log output
Date_time <- Sys.time()
log_data <- as.data.frame(Date_time, row.names = NULL)
log_data$Project_name <- project_name
log_data$Samples <- do.call(paste, c(as.list(samples_v), sep = ", "))
log_data$Organism <- organism
log_data$Approach <- approach
log_data$N_cells_raw <- do.call(paste, c(as.list(as.numeric(table(Idents(rawData_combined)))), sep = ", "))
log_data$N_cells_filtered <- do.call(paste, c(as.list(as.numeric( table(Idents(data_combined)))), sep = ", "))
log_data$Filtering_parameters <- paste("nFeature_RNA:", min_nFeature, "& log10GenesPerUMI:", min_log10GenesPerUMI, "& percent_mt:", min_percent_mt, sep = " ")
log_data$PC_number_doublet <- do.call(paste, c(as.list(dim_v), sep = ", "))

write.table(log_data, file = "1_scPipeline_QC_13052022_LOG.log", sep = "\t", row.names = FALSE)


