# This is exploring the filtered data from the Atlas of Griffiths 2017 from D6.5 to D8.5
# Paper Pijuan-Sala Griffiths Göttgens et al. Nature 2019 10.1038/s41586-019-0933-9

library(Seurat)
library(cowplot)
library(dplyr)
library(Matrix)
library(ggplot2)
library(harmony)
library(scmap)
library(scater)
library(plotly)

atlas.directory <- "D:/scRNAseq_Atlases/Atlas_E-MTAB-6967.processed.1"
output.directory <- "C:/Users/nicki/Documents/PostDoc_LSCB/20-05-13_Cardiac_Gastruloids/github/scmap_pijuansala_classifier"
de_novo_marker_discovery <- FALSE

##### Open raw data #####
setwd(dir = atlas.directory)
#atlas.data <- readLines("./atlas/raw_counts.mtx",n = 5)
atlas.metadata <- read.csv("./atlas/meta.csv")
atlas.genes <- read.csv("./atlas/genes.tsv",sep="\t",header = FALSE,as.is = TRUE)
atlas.data <- readMM("./atlas/raw_counts.mtx")
atlas.data@Dimnames[[2]] <- atlas.metadata$cell
atlas.data@Dimnames[[1]] <- atlas.genes$V2

##### Create Seurat object #####
atlas <- CreateSeuratObject(counts = atlas.data, project = "Atlas", min.cells=3, min.features = 200)
atlas <- NormalizeData(atlas, normalization.method = "LogNormalize", scale.factor = 10000)
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$stage,col.name = "stage")
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$celltype,col.name = "celltype")
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$sequencing.batch,col.name = "Replicate")

# Find which genes to use for the scmap classifier:
Idents(atlas) <- atlas$celltype
#Takes overnight to compute, so store the result in a file to avoid recomputing it too often:
if(de_novo_marker_discovery){
  atlas <- FindVariableFeatures(atlas, selection.method = "vst", nfeatures = 6000)
  markers <- FindAllMarkers(object = atlas, return.thresh = 1e-5,only.pos = T,features = VariableFeatures(atlas))
  setwd(output.directory)
  write.table(x = markers, file = "GriffithAtlas_AllMarkers_forScMapclassifier.tsv",sep = "\t")
  markers.sign <- markers[markers$avg_logFC>0.5 & markers$p_val_adj<1e-5,]
}else{
  setwd(output.directory)
  markers <- read.table("GriffithAtlas_AllMarkers_forScMapclassifier.tsv",sep = "\t")
  markers.sign <- markers[markers$avg_logFC>0.5 & markers$p_val_adj<1e-5,]
}

##### Construct the scmap classifier #####
set.seed(1234567)
scmap_classifier <- vector("list", length(unique(atlas$stage)))
names(scmap_classifier) <- unique(atlas$stage)
for(stage in unique(atlas$stage)){
  Idents(atlas) <- atlas$stage
  tmp <- subset(x = atlas, idents = stage)
  Idents(tmp) <- tmp$celltype
  tmp <- subset(x = tmp, cells = colnames(tmp)[is.na(Idents(tmp))], invert = T)
  tmp <- FindVariableFeatures(object = tmp, selection.method = 'vst', nfeatures = 200)
  atlas.sce <- as.SingleCellExperiment(x = tmp)
  rowData(atlas.sce)$feature_symbol <- rownames(atlas.sce)
  atlas.sce <- atlas.sce[!duplicated(rownames(atlas.sce)), ]
  counts(atlas.sce) <- as.matrix(counts(atlas.sce))
  logcounts(atlas.sce) <- as.matrix(logcounts(atlas.sce))
  #atlas.sce <- selectFeatures(atlas.sce,suppress_plot = F,n_features = 1000) # Rather use our own chosen features (cell type markers) that the default (genes with abnormal dropout rate), as results are much more biologically meaningful.
  rowData(atlas.sce)$scmap_features <- FALSE
  rowData(atlas.sce)[markers.sign$gene,"scmap_features"] <- TRUE
  atlas.sce <- indexCluster(atlas.sce,cluster_col = "celltype")
  scmap_classifier[[stage]] <- metadata(atlas.sce)$scmap_cluster_index
  #Visualizations for setup purpose:
  #table(rowData(atlas.sce)$scmap_features)
  #head(metadata(atlas.sce)$scmap_cluster_index)
  #heatmap(as.matrix(metadata(atlas.sce)$scmap_cluster_index))
  #scmapCluster_results <- scmapCluster(projection = atlas.sce, index_list = list(scmapclassfier = metadata(atlas.sce)$scmap_cluster_index))
  #plot(getSankey(colData(atlas.sce)$celltype, scmapCluster_results$scmap_cluster_labs[,"scmapclassfier"],plot_height = 1200))
}

saveRDS(object = scmap_classifier, file = "scmap_classifier_1000markers.rds")

##### Facultative: Exploration of atlas data ######
# Here one can choose a subset to train a classifier only at a given stage by changing the comments, 
# or keep the whole, just removing mixed_gastrulation which is another line of study that doesn't fit with the rest (as done here):
Idents(atlas) <- atlas$stage
atlas.subset <- subset(x = atlas, idents = c("mixed_gastrulation"), invert=T)
# atlas.subset <- subset(x = atlas, idents = c("E6.5"))
# atlas.subset <- subset(x = atlas, idents = c("E6.75"))
# atlas.subset <- subset(x = atlas, idents = c("E7.0"))
# atlas.subset <- subset(x = atlas, idents = c("E7.25"))
# atlas.subset <- subset(x = atlas, idents = c("E7.5"))
# atlas.subset <- subset(x = atlas, idents = c("E7.75"))
# atlas.subset <- subset(x = atlas, idents = c("E8.0"))
# atlas.subset <- subset(x = atlas, idents = c("E8.25"))
# atlas.subset <- subset(x = atlas, idents = c("E8.5"))
# atlas.subset <- subset(x = atlas, idents = c("mixed_gastrulation"))

# Remove NA cells:
Idents(atlas.subset) <- atlas.subset$celltype
atlas.subset <- subset(x = atlas.subset, cells = colnames(atlas.subset)[is.na(Idents(atlas.subset))], invert = T)

# Visualization of this subset:
atlas.subset <- FindVariableFeatures(atlas.subset, selection.method = "vst", nfeatures = 2000)
atlas.subset <- ScaleData(atlas.subset)
atlas.subset <- RunPCA(atlas.subset,npcs = 40)
ElbowPlot(atlas.subset, ndims = 40)
dims.use <- 1:30
# To visualize without batch correction:
#atlas.subset <- RunUMAP(atlas.subset,reduction = "pca",dims=dims.use,min.dist = 1,spread = 10) # Nice for D7.0. min.dist = 2,spread = 20 for D6.5
atlas.subset <- RunHarmony(object = atlas.subset, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use, plot_convergence = T, verbose = T)
atlas.subset <- RunUMAP(atlas.subset,reduction = "harmony",dims=dims.use,min.dist = 1,spread = 4,n.neighbors = 300) # min.dist = 1,spread = 10 Nice for D7.0. min.dist = 2,spread = 20 for D6.5

# Visualize metadata as annotated by Pijuan-Sala et al.
Idents(atlas.subset) <- atlas.subset$stage
Idents(atlas.subset) <- atlas.subset$celltype
Idents(atlas.subset) <- atlas.subset$Replicate
DimPlot(object = atlas.subset,reduction = "umap", label = TRUE, pt.size = 1.5,cells = sample(Cells(atlas.subset))) # export 24-26x24

# Look at gene expression:
gene_list <- c("Sox2","Otx2","T","Eomes","Cdx2","Pou5f1","Bmp4","Amn")
FeaturePlot(object = atlas.subset, features = gene_list, pt.size = 1.5, sort.cell = TRUE)
