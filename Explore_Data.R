library(Seurat)
setwd("SETYOURDIRECTORY_TO_DATA")

# Data reload
rawdata <-Read10X(data.dir = getwd(), gene.column = 1)
rawdata[1:5,1:5]
metadata <-read.table(file = "metadata.tsv.gz", header = T, sep = "\t", row.names = 1)
head(metadata) #celltype is cluster number

# Create Seurat object, including batch-corrected PCA coordinates (i.e. harmony) and UMAP. 
SO <- CreateSeuratObject(counts = rawdata, project = "Cardiac Gastruloids", meta.data = metadata)
SO <- NormalizeData(object = SO, scale.factor = 10000)
SO[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,paste0("harmony_",1:30)]),key = "harmony_",assay = "RNA")
SO[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,c("UMAP_1","UMAP_2")]),key = "UMAP_",assay = "RNA")
SO@meta.data[,grep("harmony",colnames(SO@meta.data))] <- NULL
SO@meta.data[,grep("umap",colnames(SO@meta.data))] <- NULL

# Visualize timepoints
Idents(SO) <- SO$stage
DimPlot(SO)

# Visualize cell types
Idents(SO) <- SO$celltype
DimPlot(SO, label = T)+NoLegend()

# Plot gene expression
gene_list <- c("Myl7","Sox2","Pou5f1","T")
FeaturePlot(SO, gene_list)
