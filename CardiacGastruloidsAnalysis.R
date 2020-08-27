library(rPython)
library(dplyr)
library(matrixStats)
library(Seurat)
library(BiocGenerics)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(reticulate)
library(Rmagic)
#library(phateR)
library(scmap)
library(scater)
library(harmony)
library(plotly)
library(uwot)
#library(RColorBrewer)
library(viridis)
library(DoubletFinder)

dataset.folder <- "SET_PATH_TO_RAW_DATA"
atlas.directory <- "SET_PATH_TO_ATLAS"
analysis.folder <- "SET_PATH"
classifier.folder <- "SET_PATH/scmap_pijuansala_classifier"
Days_to_analyze <- c(4,5,6,7)
Batches_to_analyze <- c(1,2,3,4,5)

# Make some output directories in the analysis folder
setwd(dir = analysis.folder)
if(!dir.exists("Idents")){dir.create("Idents")}
if(!dir.exists("OutputTables")){dir.create("OutputTables")}
if(!dir.exists("Correlations")){dir.create("Correlations")}
if(!dir.exists("GenePlots_Aligned")){dir.create("GenePlots_Aligned")}
if(!dir.exists("GenePlots_GasAlone")){dir.create("GenePlots_GasAlone")}
if(!dir.exists("DataExportGasAlone")){dir.create("DataExportGasAlone")}
if(!dir.exists("Closeups")){dir.create("Closeups")}

##### Load custom color tables #####
colors.table <- read.table(file = "InputTables/ClusterColors.tsv", sep="\t", header = T, comment.char = "", as.is = T)
colors.use_ab.initio   <- setNames(colors.table$color[!is.na(colors.table$ab_initio_identity)],colors.table$ab_initio_identity[!is.na(colors.table$ab_initio_identity)])
colors.use_transferred <- setNames(colors.table$color[!is.na(colors.table$transferred_identity)],colors.table$transferred_identity[!is.na(colors.table$transferred_identity)])
colors.use_stages <- setNames(c('#6165fb','#f854ee','#f85454','#f8ea54','#4144d4','#6165fb','#8c61fb','#b861fb','#f854ee','#f8548a','#f85554','#f87d54','#f8ea54'),c("Day4","Day5","Day6","Day7","E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"))

##### Data loading and QC - including a first ID transfer to be able to remove doublets #####
# Cardiac gastruloids loading
setwd(dir = dataset.folder)
datasets.all <- read.table(file = "DatasetsMetadata.tsv", sep = "\t",header = TRUE, stringsAsFactors = F)
datasets.all <- datasets.all[,1:(ncol(datasets.all)-2)]
datasets.all
datasets <- datasets.all[datasets.all$Day %in% Days_to_analyze & datasets.all$batch %in% Batches_to_analyze,]
rownames(datasets) <- 1:nrow(datasets)
datasets
raw.data <- list()
for(i in 1:nrow(datasets)){
  raw.data[[i]] <- Read10X(data.dir = datasets$Filename[i])
}

##### Load the cell type classifier trained on the in vivo atlas from Pijuan-Sala Griffiths Gottgens et al. Nature 2019 10.1038/s41586-019-0933-9 #####
# Note: could skip this part and do the ID transfer directly on the harmony-aligned datasets possibly. 
# Note2: Detailed examination of marker genes shows that scmap was more accurate than harmony, so keep this transfer. 
scmap_classifier <- readRDS(file = paste(classifier.folder,"scmap_classifier_1000markers.rds",sep="/"))
ref_stages <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")

# Cardiac gastruloids Seurat object
setwd(dir = analysis.folder)
SO.list <- list()
homotypic.proportion <- list()
doublet_pct <- list()
pK <- list()
BCmetric <- list()
nDoublets <- list()
nDoublets_nonhomo <- list()
for(i in 1:nrow(datasets)){
  SO.list[[i]] <- CreateSeuratObject(counts = raw.data[[i]],project = datasets$Filename[i],assay = "RNA",min.cells = 3,min.features = 100)
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Filename[i], col.name = "Dataset")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Day[i], col.name = "Day")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$batch[i], col.name = "batch")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = colSums(SO.list[[i]][grep(pattern = "^mt-", x = rownames(SO.list[[i]]), value = TRUE), ])/colSums(SO.list[[i]])*100, col.name = "percent.mito")
  
  #ggplot(SO.list[[i]]@meta.data) + geom_point(aes(x=nFeature_RNA, y=percent.mito, color=nCount_RNA))+ scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
  SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > datasets$Min.Genes[i] & percent.mito > datasets$Min.Pt.Mito[i] & percent.mito < datasets$Max.Pt.Mito[i])
  SO.list[[i]] <- NormalizeData(object = SO.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  #SO.list[[i]] <- FindVariableFeatures(object = SO.list[[i]],assay = "RNA",selection.method = "vst",nfeatures = 2000,verbose = TRUE)

  # Mini dimensionality reduction on individual datasets just for the doublet finder.
  npcs_df <- 30 # number of principal components used for doublet finder
  SO.list[[i]] <- FindVariableFeatures(SO.list[[i]],nfeatures = 1000)
  SO.list[[i]] <- ScaleData(SO.list[[i]])
  SO.list[[i]] <- RunPCA(SO.list[[i]],npcs = npcs_df)
  SO.list[[i]] <- RunUMAP(SO.list[[i]],dims = 1:npcs_df)
  
  # Apply cell type classifier (needed already now, early on, to estimate homotypic doublet proportion)
  set.seed(1234567)
  sce <- as.SingleCellExperiment(x = SO.list[[i]])
  rowData(sce)$feature_symbol <- rownames(sce)
  #sce <- sce[!duplicated(rownames(sce)), ] # Takes 20 GB of RAM and some time (probably using non-sparse matrix intermediate...)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  scmapCluster_results <- scmapCluster(projection = sce, index_list = scmap_classifier[ref_stages],threshold = 0)
  #plot(getSankey(colData(sce)$celltype, scmapCluster_results$scmap_cluster_labs[,'atlas'],plot_height = 1200))
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = scmapCluster_results$combined_labs, col.name = "celltype")
  Idents(SO.list[[i]]) <- SO.list[[i]]$celltype
  
  # # DoubletFinder: best pK identification
  # sweep.res.list <- paramSweep_v3(SO.list[[i]], PCs = 1:npcs_df, sct = FALSE)
  # sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  # bcmvn <- find.pK(sweep.stats)
  # BCmetric[[i]] <- bcmvn$BCmetric
  # ggplot(bcmvn,mapping = aes(x=pK,y=BCmetric)) + geom_point()
  # pK[[i]] <- as.numeric(as.character(bcmvn$pK))[bcmvn$BCmetric==max(bcmvn$BCmetric)] # pK corresponding to the max in BC metric
  
  # Estimate number of doublets
  homotypic.proportion[[i]] <- modelHomotypic(SO.list[[i]]$celltype)
  doublet_pct[[i]] <- ncol(SO.list[[i]])/1000*0.76 # The estimate of 10x genomics as per the RNA profiling manual are 0.76% doublet per thousand recovered cells.
  nDoublets[[i]] <- round(doublet_pct[[i]]/100*length(colnames(SO.list[[i]])))
  nDoublets_nonhomo[[i]] <- round(nDoublets[[i]]*(1-homotypic.proportion[[i]]))
}

# Check the BCmetrics, that can be used to automate the choice of the pK:
metric_avg <- BCmetric[[1]]*0
pK_ <- as.numeric(as.character(bcmvn$pK))
for(i in 1:nrow(datasets)){
  metric_avg <- metric_avg+BCmetric[[i]]/nrow(datasets)
  plot(pK_,BCmetric[[i]],type="l",col=c("red","green","blue","cyan","magenta","grey","yellow","black")[i])
  par(new=TRUE)
}
plot(pK_,metric_avg)

## Run DoubletFinder
# First try using pK = pK[[i]] but not robust enough -> datasets 5,6,8 have no clear max and end up with a random very high pK removing a whole true single cell cluster instead of isolated doublets
# We can have rather small true clusters, so we need rather low pK. The average BCmetric curve shows three local peaks: 0.02 to 0.04, 0.15-0.17, and ~0.23-0.25. Keep the first peak, which is the correct one, so choose pK=0.03 for all. 
pK.use <- 0.03
for(i in 1:nrow(datasets)){
  SO.list[[i]]@meta.data[,grep("DF.classifications",colnames(SO.list[[i]]@meta.data),value=T)] <- NULL
  SO.list[[i]]@meta.data[,grep("pANN",colnames(SO.list[[i]]@meta.data),value=T)] <- NULL
  SO.list[[i]] <- doubletFinder_v3(SO.list[[i]], PCs = 1:npcs_df, pN = 0.25, pK = pK.use, nExp = nDoublets_nonhomo[[i]], reuse.pANN = FALSE, sct = FALSE)
  FeaturePlot(SO.list[[i]],pt.size=1.5,feature=grep("pANN",colnames(SO.list[[i]]@meta.data),value=T))
  DimPlot(SO.list[[i]],pt.size = 1.5,group.by = grep("DF.classifications",colnames(SO.list[[i]]@meta.data),value=T),cols=c("red","lightgrey"))
}
rm(sce)
gc()

##### Merge all the datasets #####
names(SO.list) <- datasets$Filename
SO <- merge(x = SO.list[[1]],y = SO.list[2:length(SO.list)],add.cell.ids = datasets$Filename, project = "Cardiac_Gastruloids", merge.data = T)
SO <- AddMetaData(object = SO, metadata = SO$batch %% 2, col.name = "Replicate")
rm(raw.data)
rm(SO.list)

# Create metadata columns with unified doublet annotations accross datasets
tmp <- character(length=ncol(SO))
names(tmp) <- colnames(SO)
for(i in grep("DF.classifications",colnames(SO@meta.data),value=T)){
  j <- !is.na(SO@meta.data[,i])
  tmp[j] <- SO@meta.data[j,i]
}
SO.align <- AddMetaData(object = SO,metadata = tmp,col.name = "doublets")
tmp <- numeric(length=ncol(SO))
names(tmp) <- colnames(SO)
for(i in grep("pANN",colnames(SO@meta.data),value=T)){
  j <- !is.na(SO@meta.data[,i])
  tmp[j] <- SO@meta.data[j,i]
}
SO <- AddMetaData(object = SO.align,metadata = tmp,col.name = "pANN")

# Remove doublets
Idents(SO) <- SO$doublets
SO <- subset(x = SO, idents = "Doublet", invert = T)

# Remove spurrious columns from SO: 
for(i in grep("DF.classifications",colnames(SO@meta.data),value=T)){
  SO@meta.data[,i] <- NULL
}
for(i in grep("pANN",colnames(SO@meta.data),value=T)){
  if(i!="pANN"){
    SO@meta.data[,i] <- NULL
  }
}
SO$doublets <- NULL

# Atlas loading
setwd(dir = atlas.directory)
atlas.metadata <- read.csv("./atlas/meta.csv")
atlas.genes <- read.csv("./atlas/genes.tsv",sep="\t",header = FALSE,as.is = TRUE)
atlas.data <- readMM("./atlas/raw_counts.mtx")
atlas.data@Dimnames[[2]] <- atlas.metadata$cell
atlas.data@Dimnames[[1]] <- atlas.genes$V2

# Atlas Seurat object
atlas <- CreateSeuratObject(counts = atlas.data, project = "Atlas", min.cells=3, min.features = 200)
atlas <- NormalizeData(atlas, normalization.method = "LogNormalize", scale.factor = 10000)
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$stage,col.name = "stage")
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$celltype,col.name = "celltype")
atlas <- AddMetaData(object = atlas,metadata = atlas.metadata$sequencing.batch,col.name = "Replicate")
atlas <- AddMetaData(object = atlas,metadata = colSums(atlas[grep(pattern = "^mt-", x = rownames(atlas), value = T), ])/colSums(atlas)*100, col.name = "percent.mito")
rm(atlas.data)

# Define a subset excluding from the atlas "mixed_gastrulation", and "NA" cells. Also exclude extraembryonic tissue, not reproduced/studied in gastruloids.
Idents(atlas) <- atlas$stage
atlas.subset <- subset(x = atlas, idents = c("mixed_gastrulation"), invert=T)
Idents(atlas.subset) <- atlas.subset$celltype
atlas.subset <- subset(x = atlas.subset, cells = colnames(atlas.subset)[is.na(Idents(atlas.subset))], invert = T)
atlas.subset <- subset(atlas.subset,idents=unique(grep("ExE",atlas.subset$celltype,value = T)), invert=T)
atlas.subset <- subset(atlas.subset,idents="Parietal endoderm", invert=T)
rm(atlas)

# Clean up the memory from the raw data and get into the output folder to export whatever analysis is done. 
setwd(dir = analysis.folder)
memory.size()
gc(verbose=T)
memory.size()

#Other option: using the union of variable genes from the individual datasets. But might miss some genes which vary depending on the stage. 
#VariableFeatures(SO) <- unique(unlist(lapply(SO[datasets$OTvsPDOvsCAF=="OT"],VariableFeatures),use.names=FALSE))

##### Homogenize the metadata between SO and the atlas.subset #####
atlas.subset <- AddMetaData(atlas.subset,metadata = paste(atlas.subset$stage,"_batch",atlas.subset$Replicate,sep=""),col.name = "Dataset")
SO <- AddMetaData(SO,metadata = paste("Day",SO$Day,sep=""),col.name = "stage")
SO$Day <- NULL
SO$model <- "Gastruloids"
atlas.subset$model <- "Embryo"
atlas.subset$batch <- atlas.subset$Replicate
atlas.subset$Replicate <- atlas.subset$Replicate+max(SO$Replicate) # Make the replicate indices different

##### Variable genes - in each dataset and in the intersection, used for integration #####
n.var.feats <- 2000
atlas.subset <- FindVariableFeatures(atlas.subset,nfeatures = n.var.feats)
LabelPoints(plot = VariableFeaturePlot(atlas.subset), points = VariableFeatures(atlas.subset), repel = FALSE, xnudge = 0, ynudge = 0,size=4)
var.feats.atlas <- VariableFeatures(atlas.subset)
SO <- FindVariableFeatures(SO,nfeatures = n.var.feats)
LabelPoints(plot = VariableFeaturePlot(SO), points = VariableFeatures(SO), repel = FALSE, xnudge = 0, ynudge = 0,size=4)
var.feats.SO <- VariableFeatures(SO)
print(paste("Shared variable features (%): ",length(intersect(var.feats.SO,var.feats.atlas))/n.var.feats*100))
var.feats.merged <- intersect(var.feats.atlas,var.feats.SO)
print(paste("Number of shared variable features used for integration: ",length(var.feats.merged)))

##### Alignment procedure #####
# Atlas-PCA in shared-variable-genes space
VariableFeatures(atlas.subset) <- var.feats.merged
atlas.subset <- ScaleData(object = atlas.subset,features = var.feats.merged,do.scale=T,do.center=T)
atlas.subset <- RunPCA(object = atlas.subset,reduction.name = "pca2",reduction.key="PC2_")
ElbowPlot(object = atlas.subset,ndims = 50,reduction = "pca2")
DimPlot(atlas.subset,reduction = "pca2",label = T)+NoLegend()
dims.use <- 1:30

# Gastruloid projection in this atlas-PCA space
VariableFeatures(SO) <- var.feats.merged
pcs <- atlas.subset@reductions$pca2@feature.loadings[,dims.use]
SO <- ScaleData(SO,features = var.feats.merged,do.scale = T,do.center = T)
tmp <- t(GetAssayData(SO,slot = "scale.data")[var.feats.merged,]) %*% pcs
SO[["pca2"]] <- CreateDimReducObject(embeddings = tmp, key = "PC2_",assay = "RNA")

# Check that indeed the data is in the same PCA space:
Idents(SO) <- SO$celltype
DimPlot(SO,reduction = "pca2",label = T)+NoLegend()
gc()

# Perform batch correction on the atlas and on the gastruloids on their own, fixing one of the datasets to stay in the "pca2" space: 
atlas.subset <- RunHarmony(object = atlas.subset,group.by.vars = "Replicate",reduction = "pca2",reduction.save = "harmony",theta = 2,lambda = 1,sigma = 0.1, nclust = 200,dims.use = dims.use,plot_convergence = T, verbose = T,reference_values = 3) # Previous value: 200
SO           <- RunHarmony(object = SO          ,group.by.vars = "Replicate",reduction = "pca2",reduction.save = "harmony",theta = 2,lambda = 1,sigma = 0.1, nclust = 200,dims.use = dims.use,plot_convergence = T, verbose = T,reference_values = 0) # Previous value: 200

# Look cell numbers by cell type in the gastruloids vs embryo
Idents(SO) <- SO$celltype
Idents(atlas.subset) <- atlas.subset$celltype
celltype.numbers.atlas <- list()
celltype.numbers.SO <- list()
for(i in sort(unique(as.character((atlas.subset$celltype))),decreasing = T)){
  celltype.numbers.atlas[[i]] <- sum(atlas.subset$celltype==i)
  celltype.numbers.SO[[i]] <- sum(SO$celltype==i)
}
celltype.numbers.atlas <- unlist(celltype.numbers.atlas)
celltype.numbers.SO <- unlist(celltype.numbers.SO)
par(las=2)
par(mar=c(5,14,4,2)) 
barplot(rbind(celltype.numbers.atlas,celltype.numbers.SO),beside = T,col = c("Red","Blue"),axes = T,horiz = T)

# Exclude the cell types that are missing in the gastruloids according to scmap before running harmony: 
min.number.of.cells.per.type <- 2
excluded.cells <- subset(x = atlas.subset, idents = names(celltype.numbers.SO[celltype.numbers.SO<min.number.of.cells.per.type]), invert=F)
atlas.subset.2 <- subset(x = atlas.subset, idents = names(celltype.numbers.SO[celltype.numbers.SO<min.number.of.cells.per.type]), invert=T)

# Merging
SO.align <- merge(x = SO, y = atlas.subset.2)
SO.align[["harmony"]] <- CreateDimReducObject(embeddings = rbind(SO[["harmony"]]@cell.embeddings[,dims.use],atlas.subset.2[["harmony"]]@cell.embeddings[,dims.use]), key = "harmony_", assay = "RNA")
SO.align <- SetAssayData(object = SO.align,assay = "RNA",slot = "scale.data",new.data = cbind(GetAssayData(SO[var.feats.merged,],slot = "scale.data"),GetAssayData(atlas.subset.2[var.feats.merged,],slot = "scale.data")))
gc()

### Need to activate the option "alignment vs an atlas", which in the first step is the embryo. 
### For the second step, can chose one of the embryo datasets as being the non-moving atlas. 
### After alignment is done, need to reinclude the other cell types from the atlas ! 

# Alignment (correct for gastruloid vs embryo, and for sequencing batches/replicates, defined in a non time-biased way):
SO.align <- RunHarmony(object = SO.align,group.by.vars = "model",reduction = "harmony",reduction.save = "harmony2",theta = 2,lambda = 1,sigma = 0.1, nclust = 500,dims.use = dims.use,plot_convergence = T, verbose = T,reference_values = "Embryo") # Analysis 5: 500

# Re-include the cells that were excluded temporarily while doing the alignment. Since the alignment kept the reference atlas fixed, harmony coordinates are still comparable to PC2. 
tmp <- SO.align[["harmony2"]]@cell.embeddings[,dims.use]
tmp2 <- GetAssayData(SO.align,slot = "scale.data")
SO.align <- merge(x=SO.align,excluded.cells)
SO.align[["harmony2"]] <- CreateDimReducObject(embeddings = rbind(tmp,excluded.cells[["harmony"]]@cell.embeddings[,dims.use]), key = "harmony2_", assay = "RNA")
SO.align <- SetAssayData(object = SO.align,assay = "RNA",slot = "scale.data",new.data = cbind(tmp2,GetAssayData(excluded.cells[var.feats.merged,],slot = "scale.data")))
rm(excluded.cells)
gc()

DimPlot(SO.align,pt.size=1,reduction = "harmony2",label=T,label.size = 8) + NoLegend()

# Do an identity transfer based on same concept as scmap (nearest cluster centroid), but this time in harmony-space to benefit from the alignment: 
# Create classifier from atlas cells
Idents(SO.align) <- SO.align$model
cells.tmp <- WhichCells(SO.align,idents = "Embryo")
tmp <- SO.align@reductions$harmony2@cell.embeddings[cells.tmp,dims.use]
idents.tmp <- SO.align$celltype[cells.tmp]
centroids <- data.frame(row.names = unique(idents.tmp))
for(i in rownames(centroids)){
  for(j in colnames(tmp)){
    centroids[i,j] <- mean(tmp[idents.tmp==i,j])
  }
}
# Apply classifier
cells.tmp <- WhichCells(SO.align,idents = "Gastruloids")
tmp <- SO.align@reductions$harmony2@cell.embeddings[cells.tmp,dims.use]
idents.tmp <- SO.align$celltype[cells.tmp]
idents.tmp[] <- ""
dist <- list()
for(i in names(idents.tmp)){
  for(j in rownames(centroids)){
    dist[[j]] <- cor(x = as.numeric(tmp[i,]),y = as.numeric(centroids[j,]),method = "pearson")
  }
  idents.tmp[i] <- names(which.max(dist))
}
SO.align$orig.celltype <- SO.align$celltype
SO.align$aligned.celltype <- SO.align$celltype
SO.align$aligned.celltype[cells.tmp] <- idents.tmp
SO.align$celltype <- SO.align$aligned.celltype
Idents(SO.align) <- SO.align$celltype
DimPlot(SO.align,cells=cells.tmp,label=T,pt.size = 1.5)+NoLegend()

# could consider pre-computing nearest neighbours to save some computation time. 
# Note there is a y option for supervised dimension reduction. 

##### Scaling and PCA of atlas alone, used to generate a reference layout #####
atlas.subset <- ScaleData(object = atlas.subset,features = var.feats.atlas,do.scale = T,do.center = T)
atlas.subset <- RunPCA(object = atlas.subset,reduction.name = "pca1",reduction.key = "PC1_")
ElbowPlot(atlas.subset,ndims=50)

# Define a layout based on a tractable small subset from the atlas, with 200 cells per cluster:
setwd(analysis.folder)
seed <- 157
spread <- 2
min_dist <- 1
set.seed(seed)
cells.use <- WhichCells(object = atlas.subset,downsample = 200) #colnames(atlas.subset) #sample(colnames(atlas.subset),1000)
nn=300
local_connectivity=1 # Tried 2 and was not so convincing. Should not be more than the local intrinsic dimension of the manifold.
fast_sgd <- F # Should set it to FALSE ultimately, to get exactly reproducible results, but TRUE can be useful as it is faster for early exploration.
umap_init <- "spectral" # "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates.
reduction.use <- "pca1"
set.seed(seed)
tmp <- umap(X = Embeddings(atlas.subset[[reduction.use]])[cells.use,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T)
tmp2 <- 0*Embeddings(atlas.subset[["pca1"]])[,1:2]
tmp2[cells.use,] <- tmp$embedding
atlas.subset[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
Idents(atlas.subset) <- atlas.subset$celltype
png(filename = paste("Idents/Layout_spread",spread,"_mind",min_dist,"_seed",seed,".png",sep=""),width = 1000,height = 1000)
print(DimPlot(atlas.subset,pt.size=1,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
dev.off()

# Pick up the average position of each cluster, as a reference layout:
Idents(atlas.subset) <- atlas.subset$celltype
x_ref <- data.frame(row.names = unique(atlas.subset$celltype))
for(i in unique(atlas.subset$celltype)){
  x_ref[i,1] <- mean(atlas.subset[["umap"]]@cell.embeddings[intersect(cells.use,names(atlas.subset$celltype[atlas.subset$celltype==i])),1])
  x_ref[i,2] <- mean(atlas.subset[["umap"]]@cell.embeddings[intersect(cells.use,names(atlas.subset$celltype[atlas.subset$celltype==i])),2])
}

# Generate a random position around the layout defined positions for all cells according to their cell types:
noise <- 5
set.seed(seed)
x_ini <- data.frame(row.names = colnames(SO.align))
x_ini[,1] <- runif(length(colnames(SO.align)))*noise
x_ini[,2] <- runif(length(colnames(SO.align)))*noise
Idents(SO.align) <- SO.align$celltype
for(i in unique(SO.align$celltype)){
  print(i)
  if(i %in% rownames(x_ref)){
    x_ini[WhichCells(SO.align,idents=i),1] <- x_ini[WhichCells(SO.align,idents=i),1]+x_ref[i,1]
    x_ini[WhichCells(SO.align,idents=i),2] <- x_ini[WhichCells(SO.align,idents=i),2]+x_ref[i,2]
  }
}

# Do an integrated (aligned) umap initialized on these layout+noise positions:
for(min_dist in c(8)){ # Use these loops to screen parameters, in order to find a good umap view
  for(spread in c(25)){
    for(nn in c(300)){
      cells.use <- colnames(SO.align) #WhichCells(object = tmp,downsample = 300) #colnames(SO) #sample(colnames(SO),1000)
      
      local_connectivity=1 # Should not be more than the local intrinsic dimension of the manifold. I would have imagined 2-3 could be reasonable, but doesn't give good results. 
      fast_sgd <- T # Should set it to false ultimately, to get exactly reproducible results, but can use T to get faster for early exploration. 
      umap_init <- as.matrix(x_ini[cells.use,]) # "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
      set.seed(seed)
      reduction.use <- "harmony2"
      tmp <- umap(X = Embeddings(SO.align[[reduction.use]])[cells.use,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T)
      tmp2 <- 0*Embeddings(SO.align[["harmony2"]])[,1:2]
      tmp2[cells.use,] <- tmp$embedding
      SO.align[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
      Idents(SO.align) <- SO.align$model
      Idents(SO.align) <- SO.align$celltype
      DimPlot(SO.align,pt.size=1,reduction = "umap",label=T,label.size = 8) + NoLegend()
      
      # Plot identities on the aligned gastruloids and embryo. For final exports: 16 inches height
      for(model_plot in c("Both","Embryo","Gastruloids")){
        for(ident_plot in c("celltype","stage","model","Dataset","Replicate")){
          png(filename = paste("Idents/Aligned_spread",spread,"_mind",min_dist,"_nn",nn,"_",ident_plot,"_",model_plot,".png",sep=""),width = 1000,height = 1000)
          Idents(SO.align) <- SO.align$model
          if(model_plot=="Both"){
            cells.plot <- cells.use
          }else{
            cells.plot <- intersect(WhichCells(object = SO.align,idents = model_plot),cells.use)
          }
          Idents(SO.align) <- factor(x = SO.align@meta.data[,ident_plot],levels = sort(unique(SO.align@meta.data[,ident_plot])),ordered = T)
          if(ident_plot=="celltype"){
            print(DimPlot(SO.align,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.plot),label.size = 8, cols = colors.use_transferred[levels(Idents(SO.align))]) + NoLegend())
          }else if(ident_plot=="stage"){
            print(DimPlot(SO.align,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot), cols = colors.use_stages[levels(Idents(SO.align))])+NoLegend())
          }else if(ident_plot=="model"){
            print(DimPlot(SO.align,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot), cols = c("#eed239","#7070fc"))) #c("#e1b736","#5555ee")
          }else{
            print(DimPlot(SO.align,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot)))
          }
          dev.off()
        }
      }
    }
  }
}

#######################################################################################
###### Projection of the gastruloids alone in their own var.genes and PCA space #######
#(but initializing UMAP on the aligned embedding to provide similarly oriented views) #
#######################################################################################
cells.use <- colnames(SO.align)[SO.align$model=="Gastruloids"]
SO <- SO.align[,cells.use]
SO[["umap2"]] <- CreateDimReducObject(embeddings = Embeddings(SO[["umap"]]), key = "UMAP2_", assay = "RNA")
SO <- FindVariableFeatures(SO, nfeatures = 1000)
SO <- ScaleData(object = SO,features = VariableFeatures(SO),do.scale = T,do.center = T)
SO <- RunPCA(SO,npcs = 50,reduction.name = "pca1",reduction.key = "PC1_")
ElbowPlot(SO,ndims=50,reduction = "pca1")
DimPlot(SO,dims = c(25,26),reduction = "pca1") # Looks quite coherent until around 25. Then still a bit of info until 30. 
dims.use <- 1:30
gc()
set.seed(1)
SO <- RunHarmony(object = SO,group.by.vars = "Replicate",reduction = "pca1",reduction.save = "harmony",theta = 2,lambda = 1,sigma = 0.1,dims.use = dims.use,plot_convergence = T, verbose = T) #, nclust = 200
# ab initio reclustering of gastruloids:
SO <- FindNeighbors(object = SO, reduction = "harmony",dims = dims.use)
res <- 4.6
seed2 <- 57301
SO <- FindClusters(object = SO, resolution = res,random.seed = seed2)
SO$louvain_clusters <- SO$seurat_clusters
SO$seurat_clusters <- NULL
SO$RNA_snn_res.4.6 <- NULL
Idents(SO) <- SO$louvain_clusters
# Initialize on aligned coordinates to keep a constant orientation, trying giving it some noise to give it some more freedom to re-converge differently. In the end, more important is to increase the n_epochs to let it reshape and re-converge. 
for(init_sdev in c(1e-4)){
  for(noise_level in c(0)){
    for(min_dist in c(6)){
      for(spread in c(15)){
        nn=300
        fast_sgd <- F # Should set it to false ultimately, to get exactly reproducible results, but faster for early exploration. 
        umap_init <- as.matrix(Embeddings(SO[["umap2"]]))  # Can be a keyword like "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
        set.seed(1)
        umap_init <- umap_init + cbind(rnorm(nrow(umap_init),sd = noise_level),rnorm(nrow(umap_init),sd = noise_level))
        set.seed(1)
        reduction.use <- "harmony" # pca1 or harmony
        tmp <- umap(X = Embeddings(SO[[reduction.use]])[cells.use,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,local_connectivity=1,ret_model=T,verbose = T,n_epochs = 1000,learning_rate = 1,init_sdev = init_sdev)
        tmp2 <- 0*Embeddings(SO[["pca1"]])[,1:2]
        tmp2[cells.use,] <- tmp$embedding
        SO[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
        png(filename = paste("Noise",noise_level,"_MinDist",min_dist,"_Spread",spread,".png",sep=""),width = 800,height = 800)
        print(DimPlot(SO,pt.size=1,reduction = "umap",label=T,label.size = 8,cells=sample(colnames(SO))) + NoLegend())
        dev.off()
      }
    }
  }
}
Idents(SO) <- SO$celltype
Idents(SO) <- SO$Replicate
Idents(SO) <- SO$louvain_clusters
DimPlot(SO,pt.size=1,reduction = "umap",label=T,label.size = 8,cells=sample(colnames(SO))) + NoLegend()

# Find markers for the Louvain clusters: 
if(F){
  Idents(SO) <- SO$louvain_clusters
  dropout.rate <- rowSums(GetAssayData(SO,slot = "data")==0)/length(colnames(SO))*100
  hist(dropout.rate,100)
  tmp <- rownames(SO)[dropout.rate>50]
  markers <- FindAllMarkers(object = SO[tmp,],min.pct = 0.1,logfc.threshold = 0.22,only.pos = F)
  markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(1.5),]
  markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
  markers.sign <- markers.sign[order(markers.sign$cluster),]
  write.table(x=markers.sign,file = "./OutputTables/Gastruloid_Markers_pP3_FC1.5.tsv",sep = "\t",row.names=F,col.names = T)
}

# Assign identities to the Louvain clusters (as inputed in the Clusters Naming table) and export all the identity plots:
tmp <- read.table(file = "InputTables/ClustersNaming.tsv",sep="\t",header = T)
Idents(SO) <- SO$louvain_clusters
levels(Idents(SO)) <- tmp$celltype
SO$celltype <- Idents(SO)
for(i in c("celltype","Replicate","louvain_clusters","stage","orig.celltype","aligned.celltype")){
  Idents(SO) <- i
  if(i=="celltype"){
    png(filename = paste("Idents/GasAlone_",i,".png",sep=""),width = 800,height = 800)
    print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,label.size = 6,cells=sample(colnames(SO)),repel = T, cols = colors.use_ab.initio[levels(Idents(SO))]) + NoLegend())
    dev.off()
  }else if(i == "aligned.celltype"){
    png(filename = paste("Idents/GasAlone_",i,".png",sep=""),width = 800,height = 800)
    print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,label.size = 6,cells=sample(colnames(SO)),repel = T, cols = colors.use_transferred[levels(Idents(SO))]) + NoLegend())
    dev.off()
  }else if(i == "stage"){
    png(filename = paste("Idents/GasAlone_",i,".png",sep=""),width = 800,height = 800)
    print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,label.size = 6,cells=sample(colnames(SO)),repel = T, cols = colors.use_stages[levels(Idents(SO))]) + NoLegend())
    dev.off()
  }else{
    png(filename = paste("Idents/GasAlone_",i,".png",sep=""),width = 800,height = 800)
    print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,label.size = 6,cells=sample(colnames(SO)),repel = T) + NoLegend())
    dev.off()
  }
}
Idents(SO) <- SO$celltype
DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,label.size = 6,cells=sample(colnames(SO)),repel = T, cols = colors.use_ab.initio[levels(Idents(SO))]) + NoLegend()

library(networkD3)
plot(getSankey(SO$celltype, SO$orig.celltype, plot_height = 1000))
plot(getSankey(SO$celltype, SO$aligned.celltype, plot_height = 1000))

#####################################
##### End of the data reloading #####
#####################################

#######################################################################################
############## Explore the gastruloids aligned to Pijuan-Sala atlas ###################
#######################################################################################

##### Plot gene lists on the aligned gastruloids/embryo #####
# Find gene families
gene_tmp <- "Bmp"
gene_list <- grep(gene_tmp,rownames(SO.align),value=T)
gene_list

# Open gene lists
conn <- file("./InputTables/GeneList.txt", "r")
gene_list <- readLines(con = conn)
close(conn)

# Make custom list
gene_list <- c("Vcam1","Kdr","Sox17")

cells.use <- colnames(SO.align)
for(gene_plot in intersect(gene_list,rownames(SO.align))){
  for(model_plot in c("Embryo","Gastruloids","Both")){
    png(filename = paste("GenePlots_Aligned/",gene_plot,"_",model_plot,".png",sep=""),width = 800,height = 800)
    Idents(SO.align) <- SO.align$model
    if(model_plot=="Both"){
      cells.plot <- cells.use
    }else{
      cells.plot <- intersect(WhichCells(object = SO.align,idents = model_plot),cells.use)
    }
    Idents(SO.align) <- SO.align$celltype
    print(FeaturePlot(SO.align,features = gene_plot,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.plot),sort.cell = T)) # for random instead of sorted
    dev.off()
  }
}

Idents(SO.align) <- SO.align$celltype
DimPlot(SO.align,pt.size=1,reduction = "umap",label=T,label.size = 8) + NoLegend()

##### Comparison of markers between pairs of clusters (e.g. visceral endoderm vs gut) #####
model.use  <- "Embryo"
ident1.use <- "Allantois"
ident2.use <- "Mesenchyme"
MinRatio <- 1.5
Idents(SO.align) <- SO.align$celltype
cells.use <- colnames(SO.align)[SO.align$model==model.use]
cells.tmp1 <- WhichCells(object = SO.align, idents = ident1.use)
cells.tmp2 <- WhichCells(object = SO.align, idents = ident2.use)
cells.tmp1 <- intersect(cells.tmp1,cells.use)
cells.tmp2 <- intersect(cells.tmp2,cells.use)
markers <- FindMarkers(SO.align, ident.1 = cells.tmp1, ident.2 = cells.tmp2, only.pos = T)
markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(MinRatio),]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
rownames(markers.sign)[1:9]
writeClipboard(rownames(markers.sign))
write.table(x=markers.sign,file = paste("./OutputTables/Markers_",model.use,"_",ident1.use,ident2.use,"_R",MinRatio,".tsv",sep = ""),sep = "\t",row.names = T,col.names = NA)

##### Correlations by cell type between the embryo and the gastruloids #####
metric.use <- "cosine" # Choose spearman, pearson or cosine
avg.before.log <- FALSE # FALSE recommended. TRUE is the Seurat default. 
max.dropout.pt <- 101 # Maximum dropout rate in one or the other of the populations compared. Set to more than 100 more for no gene selection based on dropout
use.scaled.values <- TRUE # Only if avg.before.log = FALSE
genes.use <- intersect(union(var.feats.atlas,var.feats.SO),intersect(rownames(SO),rownames(atlas.subset))) # Could also consider the union of variable genes instead of the intersection. 
cell.types  <- unique(SO.align$celltype[SO.align$model == "Gastruloids"])
cells.atlas <- SO.align$model=="Embryo"
cells.gas   <- SO.align$model=="Gastruloids"
correlations <- matrix(0,ncol = length(cell.types), nrow = length(cell.types))
rownames(correlations) <- cell.types
colnames(correlations) <- cell.types
avg.atlas.by.celltype <- data.frame(matrix(0,ncol = length(cell.types), nrow = length(genes.use)))
colnames(avg.atlas.by.celltype) <- cell.types
rownames(avg.atlas.by.celltype) <- genes.use
avg.gas.by.celltype <- avg.atlas.by.celltype
dropout.atlas.by.celltype <- avg.atlas.by.celltype
dropout.gas.by.celltype   <- avg.atlas.by.celltype
# Compute the scaled data on the genes that we use:
if(use.scaled.values){
  scaled.data.atlas <- as.matrix(SO.align@assays$RNA@data[genes.use,cells.atlas])
  scaled.data.atlas <- scaled.data.atlas-matrix(rowMeans(scaled.data.atlas),nrow = dim(scaled.data.atlas)[1],ncol = dim(scaled.data.atlas)[2])
  tmp <- rowSds(scaled.data.atlas)
  tmp[tmp==0] <- 1 # To avoid division by zero of the genes which are 100% dropout in the atlas
  scaled.data.atlas <- scaled.data.atlas/matrix(tmp,nrow = dim(scaled.data.atlas)[1],ncol = dim(scaled.data.atlas)[2])
  
  scaled.data.gas <- as.matrix(SO.align@assays$RNA@data[genes.use,cells.gas])
  scaled.data.gas <- scaled.data.gas-matrix(rowMeans(scaled.data.gas),nrow = dim(scaled.data.gas)[1],ncol = dim(scaled.data.gas)[2])
  tmp <- rowSds(scaled.data.gas)
  tmp[tmp==0] <- 1 # To avoid division by zero of the genes which are 100% dropout in the gastruloids
  scaled.data.gas <- scaled.data.gas/matrix(tmp,nrow = dim(scaled.data.gas)[1],ncol = dim(scaled.data.gas)[2])
}
# Precompute all average expressions and dropout rates by cluster
for(cell.identity in cell.types){
  print(cell.identity)
  cells.atlas.tmp <- colnames(SO.align)[SO.align$celltype==cell.identity & cells.atlas]
  cells.gas.tmp   <- colnames(SO.align)[SO.align$celltype==cell.identity & cells.gas  ]
  if(avg.before.log){
    if(used.scaled.values){
      print("error, not supporting scaling on non-normalized data average")
    }
    avg.atlas <- AverageExpression(object = SO.align[,cells.atlas.tmp], assays = "RNA", features = genes.use, return.seurat = F, slot = "data")
    avg.gas   <- AverageExpression(object = SO.align[,cells.gas.tmp  ], assays = "RNA", features = genes.use, return.seurat = F, slot = "data")
    avg.atlas <- unlist(avg.atlas[[1]])
    avg.gas   <- unlist(avg.gas[[1]])
    avg.atlas <- log1p(avg.atlas)
    avg.gas   <- log1p(avg.gas)
  }else{
    if(use.scaled.values){
      avg.atlas <- setNames(rowSums(scaled.data.atlas[genes.use,cells.atlas.tmp])/length(cells.atlas.tmp),genes.use)
      if(length(cells.gas.tmp)>1){ # A few of the transferred identities have only one cell (erythroid), which needs special handling because it crashes rowSums.
        avg.gas   <- setNames(rowSums(scaled.data.gas[genes.use,cells.gas.tmp  ])/length(cells.gas.tmp),  genes.use)
      }else{
        avg.gas   <- setNames(scaled.data.gas[genes.use,cells.gas.tmp],genes.use)
      }
    }else{
      avg.atlas <- setNames(rowSums(SO.align@assays$RNA@data[genes.use,cells.atlas.tmp])/length(cells.atlas.tmp),genes.use)
      if(length(cells.gas.tmp)>1){ # A few of the transferred identities have only one cell (erythroid), which needs special handling because it crashes rowSums.
        avg.gas   <- setNames(rowSums(SO.align@assays$RNA@data[genes.use,cells.gas.tmp  ])/length(cells.gas.tmp),  genes.use)
      }else{
        avg.gas   <- setNames(SO.align@assays$RNA@data[genes.use,cells.gas.tmp],genes.use)
      }
    }
  }
  avg.atlas.by.celltype[genes.use,cell.identity] <- avg.atlas[genes.use]
  avg.gas.by.celltype[genes.use,cell.identity]   <- avg.gas[genes.use]
  dropout.atlas.by.celltype[,cell.identity] <- rowSums(SO.align[genes.use,cells.atlas.tmp]@assays$RNA@data==0)/length(cells.atlas.tmp)*100
  dropout.gas.by.celltype[,cell.identity]   <- rowSums(SO.align[genes.use,cells.gas.tmp]@assays$RNA@data==0)  /length(cells.gas.tmp  )*100
}
# Compute the correlations/similarities
for(i in 1:length(cell.types)){
  for(j in i:length(cell.types)){
    cell.identity.atlas <- cell.types[i]
    cell.identity.gas   <- cell.types[j]
    
    genes.keep <- (dropout.atlas.by.celltype[,cell.identity.atlas]<max.dropout.pt)|(dropout.gas.by.celltype[,cell.identity.gas]<max.dropout.pt)
    
    avg.atlas <- avg.atlas.by.celltype[genes.keep,cell.identity.atlas]
    avg.gas   <- avg.gas.by.celltype[genes.keep,  cell.identity.gas  ]
    if(metric.use %in% c("pearson","spearman")){
      correlation.tmp <- cor(avg.atlas,avg.gas,method=metric.use)
    }else if(metric.use == "cosine"){
      norm.tmp = sqrt(sum(avg.atlas*avg.atlas)*sum(avg.gas*avg.gas))
      correlation.tmp <- sum(avg.atlas*avg.gas)/(norm.tmp)
    }
    print(correlation.tmp)
    correlations[cell.identity.atlas,cell.identity.gas] <- correlation.tmp
  }
}
correlations2 <- correlations+t(correlations)
for(i in 1:ncol(correlations)){
  correlations2[i,i] <- correlations2[i,i]/2
}
correlations_trunc <- correlations2
correlations_trunc[correlations_trunc<0.2] <- 0  # To reorganize by similarity, only consider correlations above 0.2 (otherwise noise within the close to zero correlations makes the heatmap disorganized)
# heatmap(x = correlations2, symm = T, Rowv = as.dendrogram(hclust(proxy::dist(correlations_trunc,method="correlation"))),col = inferno(256))
# legend(x="bottomright", legend=c("min", "ave", "max"), fill=colorRampPalette(inferno(256))(3))
heatmap(x = correlations2, symm = T, Rowv = as.dendrogram(hclust(proxy::dist(correlations_trunc,method="correlation"))),col = colorRampPalette(c("lightgrey","white", "blue"))(n = 256))
legend(x="bottomright", legend=c("min", "ave", "max"), fill=colorRampPalette(c("lightgrey","white", "blue"))(n = 3))

# Histogram of correlations
hist(correlations2,50)

# Gene-gene correlations scatter-plots:
gene_correlation_closeups <- c(25,19,27,32) # Cardiomyocytes 25, gut 19 brain 27, Erythroid 32
for(i in gene_correlation_closeups){
  for(j in gene_correlation_closeups){
    cell.identity.atlas <- cell.types[i]
    cell.identity.gas   <- cell.types[j]
    
    genes.keep <- (dropout.atlas.by.celltype[,cell.identity.atlas]<max.dropout.pt)|(dropout.gas.by.celltype[,cell.identity.gas]<max.dropout.pt)
    
    avg.atlas <- avg.atlas.by.celltype[genes.keep,cell.identity.atlas]
    avg.gas   <- avg.gas.by.celltype[genes.keep,  cell.identity.gas  ]
    if(metric.use %in% c("pearson","spearman")){
      correlation.tmp <- cor(avg.atlas,avg.gas,method=metric.use)
    }else if(metric.use == "cosine"){
      norm.tmp = sqrt(sum(avg.atlas*avg.atlas)*sum(avg.gas*avg.gas))
      correlation.tmp <- sum(avg.atlas*avg.gas)/(norm.tmp)
    }
    print(paste(cell.types[i],cell.types[j],"Correlation = ",correlation.tmp))
    correlations[cell.identity.atlas,cell.identity.gas] <- correlation.tmp
    pdf(file = paste("Correlations/GasVsEmb_",gsub(pattern = "/", replacement = "", x = cell.identity.gas),"Vs",gsub(pattern = "/", replacement = "", x = cell.identity.atlas),"_rank.pdf",sep=""))
    print(ggplot(data.frame(Embryo=rank(avg.atlas),Gastruloid=rank(avg.gas),gene=genes.use[genes.keep]),aes(x=Embryo,y=Gastruloid,label=gene))+geom_point(colour="red")+geom_text()+theme_classic())
    dev.off()
    pdf(file = paste("Correlations/GasVsEmb_",gsub(pattern = "/", replacement = "", x = cell.identity.gas),"Vs",gsub(pattern = "/", replacement = "", x = cell.identity.atlas),"_scale.pdf",sep=""))
    print(ggplot(data.frame(Embryo=avg.atlas,Gastruloid=avg.gas,gene=genes.use[genes.keep]),aes(x=Embryo,y=Gastruloid,label=gene))+geom_point(colour="red")+geom_text()+theme_classic())
    dev.off()
  }
}

# Correlations between genes rather than cells: 
metric.use <- "cosine" # spearman, pearson or cosine
genes.use <- intersect(var.feats.atlas,rownames(SO))

gene.correlation <- numeric(length(genes.use))
names(gene.correlation) <- genes.use
for(gene.tmp in genes.use){
  expr.atlas <- as.numeric(avg.atlas.by.celltype[gene.tmp,])
  expr.gas   <- as.numeric(avg.gas.by.celltype[gene.tmp,])
  if(metric.use %in% c("pearson","spearman")){
    correlation.tmp <- cor(expr.atlas,expr.gas,method=metric.use)
  }else if(metric.use == "cosine"){
    norm.tmp = sqrt(sum(expr.atlas*expr.atlas)*sum(expr.gas*expr.gas))
    correlation.tmp <- sum(expr.atlas*expr.gas)/(norm.tmp)
  }
  gene.correlation[gene.tmp] <- correlation.tmp
}
avg.gene.expression <- rowSums(atlas.subset.2[genes.use,])/length(colnames(atlas.subset.2))
atlas.subset.2 <- FindVariableFeatures(atlas.subset.2)
high.var.genes.info <- HVFInfo(atlas.subset.2)
high.var.genes.info <- high.var.genes.info[genes.use,]
high.var.genes.info$avg.expr <- avg.gene.expression[genes.use]
high.var.genes.info$log1p.mean <- log1p(high.var.genes.info$mean)
high.var.genes.info$gene.name <- rownames(high.var.genes.info)
high.var.genes.info$log1p.var.std <- log1p(high.var.genes.info$variance.standardized)
highlighted_genes <- c(high.var.genes.info$gene.name[high.var.genes.info$variance.standardized>5]) #high.var.genes.info$gene.name[gene.correlation<0.3]
#highlighted_genes <- c("Hand1","Hand2","Hcn4","Tbx1","Tbx5","Isl1","Kdr","Mesp1","Gata4","Gata6")
p1 <- ggplot(high.var.genes.info,aes(log1p.var.std,gene.correlation,label=gene.name))+geom_point()
p1 + geom_text_repel(data = high.var.genes.info[highlighted_genes,],aes(log1p.var.std,gene.correlation[highlighted_genes],label=gene.name))

#######################################################################################
################# Explore the gastruloids re-analyzed alone ###########################
#######################################################################################

# Find markers by assigned cell type: 
Idents(SO) <- SO$celltype
dropout.rate <- rowSums(GetAssayData(SO,slot = "data")==0)/length(colnames(SO))*100
hist(dropout.rate,100)
tmp <- rownames(SO)[dropout.rate>50]
markers <- FindAllMarkers(object = SO[tmp,],min.pct = 0.1,logfc.threshold = 0.22,only.pos = T)
markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(1.5),]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
markers.sign <- markers.sign[order(markers.sign$cluster),]
write.table(x=markers.sign,file = "./OutputTables/Gastruloid_CellTypeMarkers_pP3_FC1.5.tsv",sep = "\t",row.names=F,col.names = T)

##### Plot gene lists on the gastruloids alone #####
# Find genes
gene_tmp <- "Irf"
gene_list <- grep(gene_tmp,rownames(SO.align),value=T)
gene_list

# Open gene lists
conn <- file("./InputTables/GeneList.txt", "r")
gene_list <- readLines(con = conn)
close(conn)

# List from all significant markers: 
gene_list <- unique(markers.sign$gene)

# Make custom list
gene_list <- c("Pax3","Pax7","Pax6","Dbx1","Dbx2")

# Make the plots
for(gene_plot in intersect(gene_list,rownames(SO))){
  png(filename = paste("./GenePlots_GasAlone/",gene_plot,".png",sep=""),width = 800,height = 800)
  print(FeaturePlot(SO,features = gene_plot,pt.size=1.5,reduction = "umap",label=T,sort.cell = T)) # ,cells = sample(cells.plot) for random instead of sorted
  dev.off()
}

# Highlight clusters
SO$highlight  <- SO$celltype=="Cardiomyocytes"
FeaturePlot(SO,"highlight")
SO$highlight2 <- SO$celltype=="Gut"
FeaturePlot(SO,c("highlight","highlight2"),blend = T)

##### Write down the processed dataset for sharing the results (DataExportGasAlone) ##### 
tmp <- GetAssayData(object = SO,slot = "counts",assay = "RNA")
write(colnames(tmp), file = "./DataExportGasAlone/barcodes.tsv")
write(rownames(tmp), file = "./DataExportGasAlone/features.tsv")
writeMM(obj = tmp, file = "./DataExportGasAlone/matrix.mtx")
tmp <- cbind(SO@meta.data, Embeddings(SO[["umap"]]), Embeddings(SO[["harmony"]])[,dims.use])
tmp$doublet_score <- tmp$pANN
tmp$pANN <- NULL
tmp$orig.ident <- NULL
tmp$orig.celltype <- NULL
tmp$pijuansala_transfer <- tmp$aligned.celltype
tmp$aligned.celltype <- NULL
tmp$highlight <- NULL
tmp$model <- NULL
tmp %>% head
write.table(x = tmp, file = "./DataExportGasAlone/metadata.tsv", sep = "\t", row.names = T, col.names = NA)

##### Export the metadata and harmony and umap coordinates for import in python / scanpy, on all cells or a close-up #####
# Name of the close-up / re-export and subset choice. Choose the first or second set or define another based on a list of louvain clusters. 
# All cells
name.tmp <- "VeloAll"
clusters_keep <- levels(SO$louvain_clusters)
is.closeup <- FALSE
# Close-up on heart development
name.tmp <- "GasHeartcloseup"
clusters_keep <- c(53,44,12,33,30,9,1,27,16,6,4,0,2,51,43,7,42,41,45,22,38,10,15,22,28)
is.closeup <- TRUE

# Computations
SO$highlight <- FALSE; SO$highlight[SO$louvain_clusters %in% clusters_keep] <- TRUE
png(filename = paste("Closeups/",name.tmp,"_CellsKept.png",sep=""),width = 800,height = 800)
FeaturePlot(SO,"highlight",sort.cell = F)
dev.off()
SO.tmp <- SO[,SO$louvain_clusters %in% clusters_keep]
if(is.closeup){ # Run only if the aim is to redo the projection rather than keep the overal view
  reduction.use <- "umap"
  SO.tmp <- FindVariableFeatures(SO.tmp)
  SO.tmp <- ScaleData(SO.tmp,features = VariableFeatures(SO.tmp))
  SO.tmp <- RunPCA(SO.tmp)
  ElbowPlot(SO.tmp)
  DimPlot(SO.tmp,reduction = "pca",dims = c(11,12)) 
  dims.use.closeup <- 1:12 # 1:12
  SO.tmp <- RunHarmony(SO.tmp,group.by.vars = "Replicate",dims.use = dims.use.closeup,nclust = 20)
  nn=50
  fast_sgd <- F # Should set it to false ultimately, to get exactly reproducible results, but faster for early exploration. 
  umap_init <- as.matrix(Embeddings(SO[,colnames(SO.tmp)][["umap"]]))  # Can be a keyword like "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
  set.seed(1)
  init_sdev <- 1e-4
  min.dist <- 6
  spread <- 12
  tmp <- umap(X = Embeddings(SO.tmp[["harmony"]])[,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,local_connectivity=1,ret_model=T,verbose = T,n_epochs = 1000,learning_rate = 1,init_sdev = init_sdev)
  tmp2 <- 0*Embeddings(SO.tmp[["pca"]])[,1:2]
  tmp2[,] <- tmp$embedding
  SO.tmp[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
  png(filename = paste(name.tmp,".png",sep=""),width = 800,height = 800)
  Idents(SO.tmp) <- SO.tmp$celltype
  print(DimPlot(SO.tmp,pt.size=2,reduction = "umap",label=T,label.size = 8,cells=sample(colnames(SO.tmp)), cols = colors.use_ab.initio[levels(Idents(SO.tmp))]) + NoLegend())
  dev.off()
  Idents(SO.tmp) <- SO.tmp$Replicate
  Idents(SO.tmp) <- SO.tmp$stage
  Idents(SO.tmp) <- SO.tmp$celltype
  Idents(SO.tmp) <- SO.tmp$louvain_clusters
  print(DimPlot(SO.tmp,pt.size = 1,label = T,reduction = reduction.use)+NoLegend())
}
getwd()
head(SO.tmp@meta.data)
write.table(x = SO.tmp@meta.data,file = paste("./OutputTables/",name.tmp,"_metadata.tsv",sep=""),append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
head(Embeddings(SO.tmp[[reduction.use]]))
write.table(x = cbind(Embeddings(SO.tmp[[reduction.use]]),Embeddings(SO.tmp[["harmony"]])), file = paste("./OutputTables/",name.tmp,"_DR.tsv",sep=""),append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
write.table(x = VariableFeatures(SO.tmp),file = paste("./OutputTables/",name.tmp,"_varfeats.tsv",sep=""),append = F,quote = F,sep = "\n",row.names = F,col.names = F)

##### Highlights on the heart close-up #####
# Mesp1 restricted to Days
Mesp1_tmp <- as.vector(GetAssayData(SO.tmp["Mesp1",]))
SO.tmp$Mesp1_Day4 <- Mesp1_tmp*(SO.tmp$stage=="Day4")
SO.tmp$Mesp1_Day5 <- Mesp1_tmp*(SO.tmp$stage=="Day5")
SO.tmp$Mesp1_Day6 <- Mesp1_tmp*(SO.tmp$stage=="Day6")
FeaturePlot(SO.tmp,c("Mesp1_Day4","Mesp1_Day5"),blend = T,pt.size = 2,sort=T)
FeaturePlot(SO.tmp,c("Mesp1_Day4","Mesp1_Day5","Mesp1_Day6","Mesp1"),min.cutoff = c(0,0,0,0),max.cutoff = c(3,3,3,3),blend = F,pt.size = 2,sort=T)

# Main Heart field signatures
SHF <- c("Tbx1", "Six2", "Six1")
FHF <- c("Tbx5","Hcn4","Nkx2-5","Hand1")
HF_markers <- c("Tbx5","Hcn4","Nkx2-5","Hand1","Tbx20","Tnnt2","Myl7","Tbx1","Isl1","Fgf8","Six1","Fgf10","Six2","Foxc1","Crabp1","Meox1","Aldh1a2","Foxf1","Sema3c","Foxc2","Prdm1","Eya1","Sall1","Tbx20")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = colSums(SO.tmp@assays$RNA@data[SHF,])/length(FHF), col.name = "SHF")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = colSums(SO.tmp@assays$RNA@data[FHF,])/length(SHF), col.name = "FHF")
FeaturePlot(object = SO.tmp, c("FHF","SHF"),blend = T,pt.size = 2,order = T,min.cutoff = 'q01',max.cutoff = 'q99')
FeaturePlot(object = SO.tmp, c("FHF","SHF","Isl1"),pt.size = 2,order = T,min.cutoff = 'q01',max.cutoff = 'q99')
FeaturePlot(object = SO.tmp, HF_markers[1:12],pt.size = 2,order = T,min.cutoff = 'q00',max.cutoff = 'q98',ncol = 4) # Export 16x24 in
FeaturePlot(object = SO.tmp, HF_markers[13:24],pt.size = 2,order = T,min.cutoff = 'q00',max.cutoff = 'q98',ncol = 4)

# Signatures for first and second Mesp1 wave: (note I had to rename the following genes to fit the current naming conventions:  "Odz4","Whsc1l1","Mll3","Fast","Flk1","Rgnef","Leprel1")
mesp_signatures <- read.table(file = "./InputTables/MespWavesSignatures_Lescroart2014.tsv", header = T, sep = "\t", as.is = T)
FW <- mesp_signatures[mesp_signatures$E6.5/mesp_signatures$E7.5>=3,] # First wave
SW <- mesp_signatures[mesp_signatures$E7.5/mesp_signatures$E6.5>=3,] # Second wave
BW <- mesp_signatures[,] # Both waves
SO.tmp <- AddMetaData(object = SO.tmp, metadata = colSums(SO.tmp@assays$RNA@data[FW$gene,]), col.name = "FW")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = colSums(SO.tmp@assays$RNA@data[SW$gene,]), col.name = "SW")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = colSums(SO.tmp@assays$RNA@data[BW$gene,]), col.name = "BW")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = (1+SO.tmp$FW)/(1+SO.tmp$SW), col.name = "WR")
FeaturePlot(object = SO.tmp, c("FW","SW","BW","FHF","SHF","Mesp1"),ncol = 3,blend = F,pt.size = 2,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))

# Signatures restricted to cells expressing Mesp1: 
Mesp1_tmp <- as.vector(GetAssayData(SO.tmp["Mesp1",]))
SO.tmp$FW_Mesp1 <- SO.tmp$FW*(Mesp1_tmp>0)
SO.tmp$SW_Mesp1 <- SO.tmp$SW*(Mesp1_tmp>0)
FeaturePlot(object = SO.tmp, c("FW","SW","BW","FW_Mesp1","SW_Mesp1","Mesp1","FHF","SHF","WR"),ncol = 3,blend = F,pt.size = 2,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))
FeaturePlot(object = SO.tmp, c("FW_Mesp1","SW_Mesp1"),ncol = 3,blend = T,pt.size = 2,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))
FeaturePlot(object = SO.tmp, c("FW_Mesp1"),blend = F,pt.size = 1,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))+NoLegend()
FeaturePlot(object = SO.tmp, c("SW_Mesp1"),blend = F,pt.size = 1,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))+NoLegend()

# Overlay of FHF and SHF signatures: 
SO.tmp$background <- 0
FeaturePlot(object = SO.tmp, c("background"),ncol = 1,blend = F,pt.size = 1,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#e5e4e4","blue"))+NoLegend()
FeaturePlot(object = SO.tmp, c("SHF"),ncol = 1,blend = F,pt.size = 1,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#FFFFFF","blue"))+NoLegend()
FeaturePlot(object = SO.tmp, c("FHF"),ncol = 1,blend = F,pt.size = 1,order = T,min.cutoff = 'q50',max.cutoff = 'q99',cols = c("#FFFFFF","red"))+NoLegend()

############################################################
##### Close-up on heart trajectories in the atlas alone ####
############################################################
# Name of the close-up
name.tmp <- "Atlas_heartcloseup"
# Subset choice
clusters_keep_atlas <- c("Intermediate mesoderm", "Mixed mesoderm", "Mesenchyme", "Pharyngeal mesoderm", "Cardiomyocytes", "Paraxial mesoderm", "Somitic mesoderm", "Epiblast", "Primitive Streak", "Caudal epiblast", "Nascent mesoderm", "Caudal Mesoderm", "Haematoendothelial progenitors", "Endothelium", "Blood progenitors 1")  # heart close-up 1
#clusters_keep <- c(0:(length(SO$louvain_clusters)-1))
# Computations
atlas.subset.3 <- SO.align[,SO.align$orig.celltype %in% clusters_keep_atlas & SO.align$model == "Embryo"]
Idents(atlas.subset.3) <- atlas.subset.3$orig.celltype
reduction.use <- "umap"
atlas.subset.3 <- FindVariableFeatures(atlas.subset.3)
atlas.subset.3 <- ScaleData(atlas.subset.3,features = VariableFeatures(atlas.subset.3))
atlas.subset.3 <- RunPCA(atlas.subset.3)
ElbowPlot(atlas.subset.3)
DimPlot(atlas.subset.3,reduction = "pca",dims = c(8,12)) 
dims.use.closeup <- 1:17 # 1:17
atlas.subset.3 <- RunHarmony(atlas.subset.3,group.by.vars = "Replicate",dims.use = dims.use.closeup,nclust = 20)
nn=300
fast_sgd <- F # Should set it to FALSE ultimately, to get exactly reproducible results, but faster with TRUE for early exploration. 
umap_init <- as.matrix(Embeddings(SO.align[,colnames(atlas.subset.3)][["umap"]]))  # Can be a keyword like "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
set.seed(1)
init_sdev <- 1e-4
min.dist <- 1
spread <- 10
tmp <- umap(X = Embeddings(atlas.subset.3[["harmony"]])[,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,local_connectivity=1,ret_model=T,verbose = T,n_epochs = 1000,learning_rate = 1,init_sdev = init_sdev)
tmp2 <- 0*Embeddings(atlas.subset.3[["pca"]])[,1:2]
tmp2[,] <- tmp$embedding
atlas.subset.3[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
png(filename = paste("Closeups/"name.tmp,".png",sep=""),width = 800,height = 800)
Idents(atlas.subset.3) <- atlas.subset.3$celltype
print(DimPlot(atlas.subset.3,pt.size=2,reduction = "umap",label=T,label.size = 8,cells=sample(colnames(atlas.subset.3)), cols = colors.use_transferred[levels(Idents(atlas.subset.3))]) + NoLegend())
dev.off()
Idents(atlas.subset.3) <- atlas.subset.3$Replicate
Idents(atlas.subset.3) <- atlas.subset.3$stage
Idents(atlas.subset.3) <- atlas.subset.3$orig.celltype
DimPlot(atlas.subset.3,pt.size = 1,label = T,reduction = reduction.use)+NoLegend()

##### 3D UMAP #####
# Compute the projection
SO <- RunUMAP(object = SO, reduction = "harmony", reduction.name = "umap3d",reduction.key = "UMAP3D_",dims = dims.use, umap.method = "uwot", min.dist = 1, n.neighbors = 300, spread = 10, seed.use = 750,n.components = 3L)

# Plot discrete properties (e.g.Dataset, celltype, Day, batch, Replicate... )
visualize <- "louvain_clusters"
plot.data <- FetchData(object = SO, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3","louvain_clusters"))
plot.data <- plot.data[sample(1:nrow(plot.data),size = nrow(plot.data)),]
colnames(plot.data) <- c("UMAP3D_1","UMAP3D_2","UMAP3D_3","louvain_clusters")
plot_ly(data = plot.data, x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, color = ~louvain_clusters, colors = rainbow(length(unique(plot.data$louvain_clusters))),type = "scatter3d", mode = "markers", marker = list(size = 2), text=~louvain_clusters,hoverinfo="text")

# Plot continuous properties (e.g. gene expression value, nFeature_RNA etc)
visualize <- "Noto" 
plot.data <- FetchData(object = SO, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3", "celltype",visualize))
plot.data <- plot.data[sample(1:nrow(plot.data),size = nrow(plot.data)),]
colnames(plot.data) <- c("UMAP3D_1","UMAP3D_2","UMAP3D_3","celltype","value")
plot_ly(data = plot.data, x = ~UMAP3D_X, y = ~UMAP3D_Y, z = ~UMAP3D_Z, color = ~value, colors = c("lightgrey","blue"),type = "scatter3d", mode = "markers", marker = list(size = 2), text=~celltype,hoverinfo="text")

##### 4D UMAP #####
SO <- RunUMAP(object = SO, reduction = "harmony", reduction.name = "umap4d",reduction.key = "UMAP4D_",dims = dims.use, umap.method = "uwot", min.dist = 0.8, spread = 15, n.neighbors = 300, seed.use = 750,n.components = 4L)
tmp <- Embeddings(SO[["umap4d"]])
tmp1 <- tmp[,1:2]
tmp2 <- tmp[,3:4]
rotation <- function(xy,theta=pi/2){
  ct <- cos(theta)
  st <- sin(theta)
  rotmat <- matrix(nrow = 2,ncol = 2,data = c(ct,st,-st,ct))
  return(t(rotmat%*%t(xy)))
}
nframes <- 1200
frames <- 1:nframes-1
getwd()
if(!dir.exists("Rotation4D")){dir.create("Rotation4D")}
for(frame in frames){
  t <- frame/nframes
  theta1 <- t*2*pi
  theta2 <- t*4*pi
  xy_proj <- data.frame(x=rotation(tmp1,theta1)[,1],y=rotation(tmp2,theta2)[,1],color=SO$celltype)
  xy_mids <- as.data.frame(matrix(data = 0,nrow = length(unique(xy_proj$color)),ncol = 2,dimnames = list(unique(xy_proj$color))))
  colnames(xy_mids) <- c("x","y")
  for(i in unique(xy_proj$color)){
    xy_mids[i,] <- colMeans(xy_proj[xy_proj$color==i,c("x","y")])
  }
  xy_mids$label <- rownames(xy_mids)
  png(filename = paste("Rotation4D/Frame_",frame,".png",sep=""),width = 1000,height = 1000)
  print(ggplot(xy_proj,aes(x = x,y = y,color=color,label=color))+geom_point()+geom_text(data=xy_mids,aes(x=x,y=y,label=label),colour="black",size=8,inherit.aes = F)+theme_light()+NoLegend())
  dev.off()
  #print(paste(frame,"_",theta1,"_",theta2,sep=""))
}

##### Data visualization and exploration #####
Idents(SO) <- SO$Dataset
Idents(SO) <- SO$celltype
Idents(SO) <- SO$Day
Idents(SO) <- SO$batch
Idents(SO) <- SO$Replicate
cells.use <- sample(colnames(SO.align)[SO.align$model=="Gastruloids"])
cells.use <- sample(colnames(SO.align)[SO.align$model=="Embryo"])
DimPlot(SO,pt.size=1.3,reduction = "umap",label=T,cells=cells.use,combine = T) # Export 16x14 to 20x14
gene_list <- c("nFeature_RNA","nCount_RNA","percent.mito","Top2a")
gene_list <- c("Myl7","T","Cd34","Noto","Pou5f1","Hoxb9","Nkx1-2","Sox17","Pdgfra","Sox2","Tubb3","Six3")
FeaturePlot(SO.align,gene_list,ncol = 4,label = T,pt.size = 2,cells = cells.use,order = T,min.cutoff = "q0",max.cutoff = "q99",reduction = "umap") # save 24x16

# Doublet scores: 
FeaturePlot(object = SO, pt.size = 1.5,features = "pANN")
Idents(SO.align) <- SO.align$model
FeaturePlot(object = SO.align, pt.size = 1.5,features = "pANN",cells=WhichCells(SO.align,idents="Gastruloids"))

# Find all markers (Embryo, Gastruloids, or both) without prefiltering on dropouts (not recommended though)
cells.use <- colnames(SO.align)
cells.use <- colnames(SO.align)[SO.align$model=="Embryo"]
cells.use <- colnames(SO.align)[SO.align$model=="Gastruloids"]
markers <- FindAllMarkers(object = SO.align[,cells.use],min.pct = 0.1,logfc.threshold = 0.25,only.pos = T)
markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(1.25),]
markers.sign[order(markers.sign$avg_logFC),]
write.table(x=markers.sign,file = "./Tables/celltypemarkersSign_embryo.tsv",sep = "\t",row.names=F,col.names = T)

##### Heatmap from canonical celltype markers on ab initio identities #####
max.cells.per.cluster <- 500
cluster.order <- read.table(file = "./InputTables/ClustersOrderHeatmap.txt",sep = "\n",header=F,as.is = T)$V1
markers <- read.table(file = "./InputTables/HighlightedMarkers.txt",sep = "\n",header=F,as.is = T)$V1
Idents(SO) <- SO$celltype
Idents(SO) <- factor(Idents(SO),levels = cluster.order)
cells.plot <- c()
for(cluster in cluster.order){
  cells.tmp <- WhichCells(object = SO,idents = cluster)
  cells.plot <- c(cells.plot,sample(cells.tmp,size = min(max.cells.per.cluster,length(cells.tmp))))
}
if(!dir.exists("Heatmaps")){dir.create("Heatmaps")}
#png(filename = paste("./Heatmaps/GasOnly_HighlightedMarkers.png",sep=""),width = 1000,height = 1000)
pdf(file = paste("./Heatmaps/GasOnly_HighlightedMarkers.pdf",sep=""), width = 16, height = 10)
print(DoHeatmap(object = SO, cells = cells.plot, features = markers, slot = "data",raster = T,draw.lines = T, lines.width = 15,group.colors = colors.use_ab.initio[levels(Idents(SO))])+ scale_fill_gradientn(colors = c("white","blue","black")))+NoLegend()
dev.off()

# Manually pick up cells to find markers
Idents(SO.align) <- SO.align$celltype
cells.use <- colnames(SO.align)[SO.align$model=="Embryo"]
plot1 <- DimPlot(object = SO.align,cells = cells.use, pt.size = 1,reduction = "umap")
cells.tmp <- CellSelector(plot1)
markers <- FindMarkers(object = SO.align, ident.1 = cells.tmp, ident.2 = setdiff(x = cells.use,y = cells.tmp))
markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(1.5),]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
head(markers.sign,50)
writeClipboard(rownames(markers.sign))
gene_list <- rownames(markers.sign)
write.table(x=markers.sign,file = "./Tables/celltypemarkers_gas.tsv",sep = "\t",row.names = F,col.names = T)

# Plot to highlight the cells selected in the QC space (to check no artefacts)
tmp <- SO.align[,cells.use]@meta.data
tmp$cells.tmp <- "Other_Cells"
tmp[cells.tmp,"cells.tmp"] <- "Highlighted_Cells"
ggplot(tmp) + geom_point(aes(x=nFeature_RNA, y=percent.mito, color=cells.tmp))+scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+scale_color_manual(values=c("Red", "Lightgrey"))+theme_light()
ggplot(tmp) + geom_violin(aes(x=cells.tmp,y=nFeature_RNA,color=cells.tmp))+ggplot2::scale_colour_manual(values=c("Red","Darkgrey"))+theme_light()

##### Aligned data export #####
# # Export the combined datasets as a loom file:
# library(reticulate)
# seuratobject_ad <- Convert(from=SO.align, to="anndata", filename="SO.align_v4.h5ad")

# Write down tables with the metadata and harmony2 and umap coordinates of the aligned datasets:
getwd()
head(SO.align@meta.data)
write.table(x = SO.align@meta.data,file = "./OutputTables/SO.align_metadata.tsv",append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
head(Embeddings(SO.align[["umap"]]))
write.table(x = cbind(Embeddings(SO.align[["umap"]]),Embeddings(SO.align[["harmony2"]])), file = "./OutputTables/SO.align_DR.tsv",append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
write.table(x = var.feats.merged,file = "./OutputTables/SO.align_var.feats.merged.tsv",append = F,quote = F,sep = "\n",row.names = F,col.names = F)

sessionInfo()
