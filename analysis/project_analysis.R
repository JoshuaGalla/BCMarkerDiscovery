library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(glue)

# Step 1: Load Data ---------------------------------

#produces a sparse matrix indicating gene expression per cell per patient:
#column = cell specific to certain patient ID (Ex: CID3586)
#row = RNA molecule count for a specific gene (i.e feature)
data <- Read10X("./GSE176078_Wu_etal_2021_BRCA_scRNASeq", gene.column = 1) #takes a couple mins
#head(data)

#nCount_RNA = Total # of RNA molecules per cell; higher number = more RNA content or more sequencing depth
#nFeature_RNA = # of unique genes expressed per cell; higher number = more diverse cellular expression
metadata <- read.csv("./GSE176078_Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
#head(metadata)

# Step 2: Setting up Seurat Object --------------------------

#Creating Seurat Object of sparse/count matrix using data and metadata - takes a couple mins
data <- CreateSeuratObject(data, project = "545_group_project", min.cells = 3, min.features = 200)

# Step 3: QC and subsetting data into HER2+, ER+, and TNBC subtypes -------

#set QC thresholds to filter out low quality reads - can be adjusted as needed
nfeature_min <- 200
nfeature_max <- 2500
ncount_min <- 500
percent_mt_max <- 10 #samples that have greater than 10-15% percent.mt represent high cellular contamination from low quality or cell death

#Adding "percent.mt" column from metadata to data object
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
#head(data@meta.data)

#add subset data type column from metadata to data
data$subtype <- metadata$subtype

#subset data based on QC counts (features/counts/percent.mt)
data <- subset(data, subset = nFeature_RNA > nfeature_min & nFeature_RNA < nfeature_max & nCount_RNA > ncount_min & percent.mt < percent_mt_max)
#head(data)

#subset data based on subtype
subtype <- "HER2+" #options include "HER2+", "ER+", or "TNBC"

data <- subset(data, subset = subtype == subtype)
#print(length(data$nCount_RNA)) #total number of cells in subset

# Step 4: Normalization: CPM ------------------

#Raw Count normalization (CPM) = (count on features)/(library size)*1,000,000
data <- NormalizeData(data, normalization.method = "RC", scale.factor = 1e6)
#head(data[["RNA"]]$counts, 20)

# Step 5: Plots ---------

#Violin plots: set features to be analyzed
markers <- c("nFeature_RNA", "nCount_RNA")

data_VlnPlot <- VlnPlot(data, features = markers, pt.size = 0.1, ncol = 2, combine = FALSE)
names(data_VlnPlot) <- markers
data_VlnPlot$nFeature_RNA & ggtitle(glue("{subtype} nFeature_RNA Raw Counts")) #nFeature_RNA count
data_VlnPlot$nCount_RNA & ggtitle(glue("{subtype} nCount_RNA Raw Counts")) #nCount_RNA counts

#Scatter plots: - title indicates R^2 value of linear fit
FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

#ER+
FeatureScatter(ER_data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

# Step 6: Highly Variable Feature (Gene Expression) Selection ---------------

#set number of highly variable features to be considered/analyzed (n_highvarfeats) and labeled (n_labels) 
n_highvarfeats = 2000
n_labels = 10
#enables overlap of labels on graph
options(ggrepel.max.overlaps = Inf) 

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = n_highvarfeats)
expr_data_high <- head(VariableFeatures(data), n_labels)
expr_data_plot <- VariableFeaturePlot(data)
LabelPoints(plot = expr_data_plot, points = expr_data_high, repel = TRUE, xnudge = 0, ynudge = 0)

# Step 7: Scaling Data (via linear transform) ---------------

data_genes <- rownames(data)
data <- ScaleData(data, features = data_genes)

# Step 8: Linear Dimension Reduction (PCA): ----------------------

#stores 5 (default) different PCAs containing largest variances (positive and negative) in gene expression in decreasing order
#ex: PCA1 contains the largest source of variation, PCA2 contains the second largest source, etc.

data <- RunPCA(data, features = VariableFeatures(object = data))

#set the number of PCA plots to be used to visualize gene expression variance
dim_min = 1
dim_max = 2

#Visualizes the top sets of genes that are associated with reduction components of PCA:
VizDimLoadings(data, dims = dim_min:dim_max, reduction = "pca")

#Graph of output of PCA where each point is a cell in position based on its reduction component:
DimPlot(data, reduction = "pca")

#Step 9: Determining dimensionality of data --------------------

#can use an elbowplot to identify how many components to include - helps remove potential technical noise
ElbowPlot(data)

#shows an "elbow" around PC 9
elbow_dims = 9 #used to determine dimensions for clustering and UMAP

#Step 10: Clustering ------------------------

#Performs graph-based clustering via K-nearest neighbor (euclidean distance)
#edges are drawn between cells of similar gene/feature expression
#clustering resolution differs for number of cells - 0.4-1.2 recommended for 3k cells:
#all subtype = 0.5 resolution (2k cells)

data <- FindNeighbors(data, dims = 1:elbow_dims) #KNN graph
data <- FindClusters(data, resolution = 0.5) #Louvain algorithm
#head(Idents(data))

#Step 11: UMAP -------------------

#preserves local distances between cell relationships and ensures cell co-localized based on gene expression
#less ideal for identifying global relationships

data <- RunUMAP(data, dims = 1:elbow_dims)
DimPlot(data, reduction = "umap", label = TRUE)

#Compares markers via differential expression to all other cells within a single cluster
data_markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#head(data_markers)

#identify largest log fold changes for upregulated genes per cluster
top_markers <- data_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

#view(top_markers)

data_markers <- FindAllMarkers(data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#head(data_markers)

#identify largest log fold changes for downregulated genes per cluster
bottom_markers <- data_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 10) %>%
  ungroup()

#view(bottom_markers)

#Step 13: Assign cell type to clusters -------------------

#Sources: Cell Marker 2.0 and The Human protein Atlas

#HER2+
HER2_ids <- c("Naive CD4+ T", "Treg", "CD8+ T Cell", "Macrophage", "Helper T", "Cytotoxic T Cell", "Stem Cell", "Chondrocyte", "?", "B Cell", "Epithelial Cell", "NK Cell", "M2 Macrophage", "Memory B cell", "Granulocyte", "Keratinocyte", "CD27+ B Cell", "Epithelial Cell")
names(HER2_ids) <- levels(data)
data <- RenameIdents(data, HER2_ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5)

#ER+
ER_ids <- c("Naive CD4+ T", "CD8+ T Cell", "Fibroblast", "Granulocyte", "Treg", "Monocyte", "Basophil", "Epithelial Cell", "B Cell", "Endothelial Cell", "T Cell", "NK CD8+ T Cell", "Mature B Cell", "B Cell", "Monocyte", "Dendritic Cell")
names(ER_ids) <- levels(data)
data <- RenameIdents(data, ER_ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5)

#TNBC
TNBC_ids <- c("Naive CD4+ T", "CD4+ T Cell", "NK Cell", "Naive B Cell", "CD8+ T Cell", "Microglia", "M2 Macrophage", "B Cell", "NK Cell", "Cancer Cell", "T Cell", "Macrophage", "Fibroblast", "Pericyte", "Stem Cell", "MK67+ Progenitor", "Epithelial Cell", "Cancer Cell", "Neural Progenitor Cell", "Naive CD4+ T", "Macrophage", "Myofibroblast", "B Cell", "Fibroblast")
names(TNBC_ids) <- levels(data)
Tdata <- RenameIdents(data, TNBC_ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5)
