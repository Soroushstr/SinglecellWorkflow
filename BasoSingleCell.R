library(Seurat)
library(dplyr)
library(Matrix)
##### sandbox #####
pbmc.data <- Read10X(data.dir = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/Lung", gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)


##### Lung HCA Naftali #####
adata.lung <- Read10X(data.dir = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/Lung", gene.column = 1)
adata <- CreateSeuratObject(counts = adata.lung, project = "lung", min.cells = 3, min.features = 200)
adata
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
head(adata@meta.data, 5)
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(adata), 10)
plot1 <- VariableFeaturePlot(adata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)
plot1 + plot2
all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
tcell <- subset(adata, subset = CD4 > 0)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
print(adata[["pca"]], dims = 1:10, nfeatures = 20)
DimHeatmap(adata, dims = 1:10, cells = 500, balanced = TRUE)
adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.5)
adata <- RunUMAP(adata, dims = 1:10)
DimPlot(adata, reduction = "umap")
adata.markers <- FindAllMarkers(adata, only.pos = TRUE)
VlnPlot(adata, features = c("FCER1A"))

adata.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
adata.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(adata, features = top10$gene)
FeaturePlot(adata, features = c("FCER1A","CD14","CCR7"))
  head(cluster2.markers, n = 5)
w  
  
tcell <- RunPCA(tcell, features = VariableFeatures(object = tcell))
print(adata[["pca"]], dims = 1:10, nfeatures = 20)
DimHeatmap(tcell, dims = 1:10, cells = 500, balanced = TRUE)
tcell <- FindNeighbors(tcell, dims = 1:10)
tcell <- FindClusters(tcell, resolution = 0.5)
tcell <- RunUMAP(tcell, dims = 1:10)
DimPlot(tcell, reduction = "umap")
tcell.markers <- FindAllMarkers(tcell, only.pos = TRUE)
VlnPlot(adata, features = c("FCER1A"))

adata.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
tcell.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(tcell, features = top10$gene)
FeaturePlot(tcell, features = c("IL1R2"))

##### COOTA Pipeline #####
adata.data = read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM3943045_Adult-Bone-Marrow1_dge.txt.gz",row.names = 1,header = T)
adata.data_500more = adata.data[,colSums(adata.data)>=300]
colnames(adata.data_500more) <- paste("2",colnames(adata.data_500more),sep = ".")
colnames(adata.data_500more) <- paste("PeriBlood",colnames(adata.data_500more),sep = "_")
adata <- CreateSeuratObject(counts = Matrix(as.matrix(adata.data_500more),sparse=T),
                            min.cells = 3, min.features = 300,names.delim = "\\.")
adata
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
head(adata@meta.data, 5)
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
adata <- subset(adata, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
DimHeatmap(adata, dims = 1:10, cells = 500, balanced = TRUE)
adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.5)
adata <- RunUMAP(adata, dims = 1:10)
DimPlot(adata, reduction = "umap")
adata.markers <- FindAllMarkers(adata, only.pos = TRUE)
VlnPlot(adata, features = c("FCER1A"))

adata.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(adata, features = top10$gene)
FeaturePlot(adata, features = c("IL7R","CD14","CD8A","GNLY","MS4A1","FCGR3A","FCER1A","PTCRA"))
FeaturePlot(adata, features = c("CST7"))
new.cluster.ids <- c("Naive CD4+ T", "CD14+ Mono", "CD8 T", "NK Cell", "B",  "FCGR3A+ Mono", "DC", "Platelet")
names(new.cluster.ids) <- levels(adata)
adata <- RenameIdents(adata, new.cluster.ids)
DimPlot(adata, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
baso <- subset(adata,subset = FCER1A >0)
baso <- RunPCA(baso, features = VariableFeatures(object = baso))
DimHeatmap(baso, dims = 1:10, cells = 500, balanced = TRUE)
baso <- FindNeighbors(baso, dims = 1:10)
baso <- FindClusters(baso, resolution = 0.5)
baso <- RunUMAP(baso, dims = 1:10)
DimPlot(baso, reduction = "umap")
baso.markers <- FindAllMarkers(baso, only.pos = TRUE)
baso.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(baso, features = top10$gene)
VlnPlot(adata, features = c("FCER1A"))
##### Large PB #####
pbmc.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008638_Adult-Peripheral-Blood1_dge.txt.gz",row.names = 1,header = T)
pbmc.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008639_Adult-Peripheral-Blood2_dge.txt.gz",row.names = 1,header = T)
pbmc.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008640_Adult-Peripheral-Blood3-1_dge.txt.gz",row.names = 1,header = T)
pbmc.3.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008641_Adult-Peripheral-Blood3-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008642_Adult-Peripheral-Blood4-1_dge.txt.gz",row.names = 1,header = T)
pbmc.4.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008643_Adult-Peripheral-Blood4-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008644_Adult-Peripheral-Blood4-3_dge.txt.gz",row.names = 1,header = T)

pbmc.1.data_500more = pbmc.1.data[,colSums(pbmc.1.data)>=300]
colnames(pbmc.1.data_500more) <- paste("1",colnames(pbmc.1.data_500more),sep = ".")
colnames(pbmc.1.data_500more) <- paste("PeriBlood",colnames(pbmc.1.data_500more),sep = "_")
pbmc.1 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.1.data_500more),sparse=T),
                            min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.1[["percent.mt"]] <- PercentageFeatureSet(pbmc.1, pattern = "^MT-")
pbmc.1 <- subset(pbmc.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.1 <- NormalizeData(pbmc.1, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.2.data_500more = pbmc.2.data[,colSums(pbmc.2.data)>=300]
colnames(pbmc.2.data_500more) <- paste("2",colnames(pbmc.2.data_500more),sep = ".")
colnames(pbmc.2.data_500more) <- paste("PeriBlood",colnames(pbmc.2.data_500more),sep = "_")
pbmc.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.2, pattern = "^MT-")
pbmc.2 <- subset(pbmc.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.2 <- NormalizeData(pbmc.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.3.data_500more = pbmc.3.data[,colSums(pbmc.3.data)>=300]
colnames(pbmc.3.data_500more) <- paste("2",colnames(pbmc.3.data_500more),sep = ".")
colnames(pbmc.3.data_500more) <- paste("PeriBlood",colnames(pbmc.3.data_500more),sep = "_")
pbmc.3 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.3.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.3[["percent.mt"]] <- PercentageFeatureSet(pbmc.3, pattern = "^MT-")
pbmc.3 <- subset(pbmc.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.3 <- NormalizeData(pbmc.3, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.3.2.data_500more = pbmc.3.2.data[,colSums(pbmc.3.2.data)>=300]
colnames(pbmc.3.2.data_500more) <- paste("2",colnames(pbmc.3.2.data_500more),sep = ".")
colnames(pbmc.3.2.data_500more) <- paste("PeriBlood",colnames(pbmc.3.2.data_500more),sep = "_")
pbmc.3.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.3.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.3.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.3.2, pattern = "^MT-")
pbmc.3.2 <- subset(pbmc.3.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.3.2 <- NormalizeData(pbmc.3.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.data_500more = pbmc.4.data[,colSums(pbmc.4.data)>=300]
colnames(pbmc.4.data_500more) <- paste("2",colnames(pbmc.4.data_500more),sep = ".")
colnames(pbmc.4.data_500more) <- paste("PeriBlood",colnames(pbmc.4.data_500more),sep = "_")
pbmc.4 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4[["percent.mt"]] <- PercentageFeatureSet(pbmc.4, pattern = "^MT-")
pbmc.4 <- subset(pbmc.4, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4 <- NormalizeData(pbmc.4, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.2.data_500more = pbmc.4.2.data[,colSums(pbmc.4.2.data)>=300]
colnames(pbmc.4.2.data_500more) <- paste("2",colnames(pbmc.4.2.data_500more),sep = ".")
colnames(pbmc.4.2.data_500more) <- paste("PeriBlood",colnames(pbmc.4.2.data_500more),sep = "_")
pbmc.4.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.4.2, pattern = "^MT-")
pbmc.4.2 <- subset(pbmc.4.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4.2 <- NormalizeData(pbmc.4.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.3.data_500more = pbmc.4.3.data[,colSums(pbmc.4.3.data)>=300]
colnames(pbmc.4.3.data_500more) <- paste("2",colnames(pbmc.4.3.data_500more),sep = ".")
colnames(pbmc.4.3.data_500more) <- paste("PeriBlood",colnames(pbmc.4.3.data_500more),sep = "_")
pbmc.4.3 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.3.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4.3[["percent.mt"]] <- PercentageFeatureSet(pbmc.4.3, pattern = "^MT-")
pbmc.4.3 <- subset(pbmc.4.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4.3 <- NormalizeData(pbmc.4.3, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.normalized <- merge(pbmc.1, y = c(pbmc.2, pbmc.3, pbmc.3.2, pbmc.4, pbmc.4.2, pbmc.4.3), 
                         add.cell.ids = c("R1", "R2", "R3", "R3.2", "R4", "R4.2", "R4.3"), project = "PBMC12K",
                         merge.data = TRUE)

pbmc.baso <- subset(pbmc.normalized,subset = IL3RA >0 & FCER1A > 0)

pbmc.normalized <- FindVariableFeatures(pbmc.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc.normalized)
pbmc.normalized <- ScaleData(pbmc.normalized, features = all.genes)
pbmc.normalized <- RunPCA(pbmc.normalized, features = VariableFeatures(object = pbmc.normalized))
DimHeatmap(pbmc.normalized, dims = 1:10, cells = 500, balanced = TRUE)
pbmc.normalized <- JackStraw(object = pbmc.normalized, reduction = "pca")
pbmc.normalized <- ScoreJackStraw(pbmc.normalized, dims = 1:20)
JackStrawPlot(object = pbmc.normalized, dims = 1:20)
ElbowPlot(pbmc.normalized)
pbmc.normalized <- FindNeighbors(pbmc.normalized, dims = 1:11)
pbmc.normalized <- FindClusters(pbmc.normalized, resolution = 0.5)
pbmc.normalized <- RunUMAP(pbmc.normalized, dims = 1:11)
DimPlot(pbmc.normalized, reduction = "umap")
FeaturePlot(pbmc.normalized, features = c("FCER1A", "IL3RA"), blend = TRUE)
pbmc.normalized2 <- JoinLayers(pbmc.normalized)
pbmc.normalized2.markers <- FindAllMarkers(pbmc.normalized2, only.pos = TRUE)

VlnPlot(pbmc.normalized2, features = c("FCER1A","IL3RA","ITGAM", "CD69"))

pbmc.normalized2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc.normalized2, features = top10$gene)

new.cluster.ids <- c("Classical Mono", "Erythroid 1", "NK Cell", "Naive T cell", "FCGR3A+ Mono", "DC",
                     "Erythroid 2", "Macrophage", "Erythroid 3", "Granulocyte", "Plasmacytoid DCs",
                     "Plasma cells", "B cells", "NKT")
names(new.cluster.ids) <- levels(pbmc.normalized)
pbmc.normalized <- RenameIdents(pbmc.normalized, new.cluster.ids)
DimPlot(pbmc.normalized, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(pbmc.normalized2, features = c("FCER1A","IL3RA","ITGAM", "CD69"))

baso <- subset(pbmc.normalized2,subset = IL3RA >0 & FCER1A > 0)
baso <- RunPCA(baso, features = VariableFeatures(object = baso))
DimHeatmap(baso, dims = 1:10, cells = 500, balanced = TRUE)
baso <- FindNeighbors(baso, dims = 1:10)
baso <- FindClusters(baso, resolution = 0.5)
baso <- RunUMAP(baso, dims = 1:10)
DimPlot(baso, reduction = "umap")
baso.markers <- FindAllMarkers(baso, only.pos = TRUE)
baso.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
VlnPlot(baso, features = c("FCER1A"))
baso.genes <- row.names(baso)
baso.genes.df <- as.data.frame(baso.genes)
FeaturePlot(baso, features = c("KIT","IL3RA","ANPEP","CCL4","CD9","IL4"))
new.cluster.ids <- c("Basophil","Plasmacytoid DC")
names(new.cluster.ids) <- levels(baso)
baso <- RenameIdents(baso, new.cluster.ids)
DimPlot(baso, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DoHeatmap(baso, features = top10$gene)

##### Large Cord Blood CD34+ #####
cb34.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008673_Cord-Blood-CD34P1_dge.txt.gz",row.names = 1,header = T)
cb34.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008671_Cord-Blood2-CD34P2-1_dge.txt.gz",row.names = 1,header = T)
cb34.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008672_Cord-Blood2-CD34P2-2_dge.txt.gz",row.names = 1,header = T)

cb34.data_500more = cb34.data[,colSums(cb34.data)>=300]
colnames(cb34.data_500more) <- paste("1",colnames(cb34.data_500more),sep = ".")
colnames(cb34.data_500more) <- paste("CordBloodCD34",colnames(cb34.data_500more),sep = "_")
cb34 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
cb34[["percent.mt"]] <- PercentageFeatureSet(cb34, pattern = "^MT-")
cb34 <- subset(cb34, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34 <- NormalizeData(cb34, normalization.method = "LogNormalize", scale.factor = 10000)

cb34.2.data_500more = cb34.2.data[,colSums(cb34.2.data)>=300]
colnames(cb34.2.data_500more) <- paste("1",colnames(cb34.2.data_500more),sep = ".")
colnames(cb34.2.data_500more) <- paste("CordBloodCD34",colnames(cb34.2.data_500more),sep = "_")
cb34.2 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
cb34.2[["percent.mt"]] <- PercentageFeatureSet(cb34.2, pattern = "^MT-")
cb34.2 <- subset(cb34.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34.2 <- NormalizeData(cb34.2, normalization.method = "LogNormalize", scale.factor = 10000)

cb34.3.data_500more = cb34.3.data[,colSums(cb34.3.data)>=300]
colnames(cb34.3.data_500more) <- paste("1",colnames(cb34.3.data_500more),sep = ".")
colnames(cb34.3.data_500more) <- paste("CordBloodCD34",colnames(cb34.3.data_500more),sep = "_")
cb34.3 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.3.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
cb34.3[["percent.mt"]] <- PercentageFeatureSet(cb34.3, pattern = "^MT-")
cb34.3 <- subset(cb34.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34.3 <- NormalizeData(cb34.3, normalization.method = "LogNormalize", scale.factor = 10000)

cd34.normalized <- merge(cb34, y = c(cb34.2, cb34.3), 
                         add.cell.ids = c("R1", "R2", "R3"), project = "CD34P",
                         merge.data = TRUE)

cd34.normalized <- FindVariableFeatures(cd34.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cd34.normalized)
cd34.normalized <- ScaleData(cd34.normalized, features = all.genes)
cd34.normalized <- RunPCA(cd34.normalized, features = VariableFeatures(object = cd34.normalized))
DimHeatmap(cd34.normalized, dims = 1:10, cells = 500, balanced = TRUE)
cd34.normalized <- JackStraw(object = cd34.normalized, reduction = "pca")
cd34.normalized <- ScoreJackStraw(cd34.normalized, dims = 1:20)
JackStrawPlot(object = cd34.normalized, dims = 1:20)
ElbowPlot(cd34.normalized)
cd34.normalized <- FindNeighbors(cd34.normalized, dims = 1:11)
cd34.normalized <- FindClusters(cd34.normalized, resolution = 0.5)
cd34.normalized <- RunUMAP(cd34.normalized, dims = 1:11)
DimPlot(cd34.normalized, reduction = "umap")
cd34.normalized2 <- JoinLayers(cd34.normalized)
cd34.normalized2.markers <- FindAllMarkers(cd34.normalized2, only.pos = TRUE)

VlnPlot(cd34.normalized2, features = c("FCER1A","IL3RA","ITGAM"))

cd34.normalized2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(cd34.normalized2, features = top10$gene)

##### Bone Marrow #####
bm.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM3943045_Adult-Bone-Marrow1_dge.txt.gz",row.names = 1,header = T)
bm.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM3980128_Adult-Bone-Marrow2_dge.txt.gz",row.names = 1,header = T)

bm.1.data_500more = bm.1.data[,colSums(bm.1.data)>=300]
colnames(bm.1.data_500more) <- paste("1",colnames(bm.1.data_500more),sep = ".")
colnames(bm.1.data_500more) <- paste("CordBloodCD34",colnames(bm.1.data_500more),sep = "_")
bm.1 <- CreateSeuratObject(counts = Matrix(as.matrix(bm.1.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
bm.1[["percent.mt"]] <- PercentageFeatureSet(bm.1, pattern = "^MT-")
bm.1 <- subset(bm.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
bm.1 <- NormalizeData(bm.1, normalization.method = "LogNormalize", scale.factor = 10000)

bm.2.data_500more = bm.2.data[,colSums(bm.2.data)>=300]
colnames(bm.2.data_500more) <- paste("1",colnames(bm.2.data_500more),sep = ".")
colnames(bm.2.data_500more) <- paste("CordBloodCD34",colnames(bm.2.data_500more),sep = "_")
bm.2 <- CreateSeuratObject(counts = Matrix(as.matrix(bm.2.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
bm.2[["percent.mt"]] <- PercentageFeatureSet(bm.2, pattern = "^MT-")
bm.2 <- subset(bm.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
bm.2 <- NormalizeData(bm.2, normalization.method = "LogNormalize", scale.factor = 10000)

bm.normalized <- merge(bm.1, y = bm.2, 
                       add.cell.ids = c("R1", "R2"), project = "BoneMarrow",
                       merge.data = TRUE)

bm.baso <- subset(bm.normalized,subset = IL3RA >0 & FCER1A > 0)

bm.normalized <- FindVariableFeatures(bm.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bm.normalized)
bm.normalized <- ScaleData(bm.normalized, features = all.genes)
bm.normalized <- RunPCA(bm.normalized, features = VariableFeatures(object = bm.normalized))
DimHeatmap(bm.normalized, dims = 1:10, cells = 500, balanced = TRUE)
bm.normalized <- JackStraw(object = bm.normalized, reduction = "pca")
bm.normalized <- ScoreJackStraw(bm.normalized, dims = 1:20)
JackStrawPlot(object = bm.normalized, dims = 1:20)
ElbowPlot(bm.normalized)
bm.normalized <- FindNeighbors(bm.normalized, dims = 1:11)
bm.normalized <- FindClusters(bm.normalized, resolution = 0.5)
bm.normalized <- RunUMAP(bm.normalized, dims = 1:11)
DimPlot(bm.normalized, reduction = "umap")
bm.normalized2 <- JoinLayers(bm.normalized)
bm.normalized2.markers <- FindAllMarkers(bm.normalized2, only.pos = TRUE)
bm.markers.df <- data.frame(bm.normalized2.markers)

FeaturePlot(bm.normalized, features = c("FCER1A", "IL3RA"), blend = TRUE)

VlnPlot(bm.normalized2, features = c("FCER1A","IL3RA"))

bm.normalized2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(bm.normalized2, features = top10$gene)

##### Spleen #####
sp.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008649_Adult-Spleen1-1_dge.txt.gz",row.names = 1,header = T)
sp.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008650_Adult-Spleen1-2_dge.txt.gz",row.names = 1,header = T)

sp.1.data_500more = sp.1.data[,colSums(sp.1.data)>=300]
colnames(sp.1.data_500more) <- paste("1",colnames(sp.1.data_500more),sep = ".")
colnames(sp.1.data_500more) <- paste("CordBloodCD34",colnames(sp.1.data_500more),sep = "_")
sp.1 <- CreateSeuratObject(counts = Matrix(as.matrix(sp.1.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
sp.1[["percent.mt"]] <- PercentageFeatureSet(sp.1, pattern = "^MT-")
sp.1 <- subset(sp.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
sp.1 <- NormalizeData(sp.1, normalization.method = "LogNormalize", scale.factor = 10000)

sp.2.data_500more = sp.2.data[,colSums(sp.2.data)>=300]
colnames(sp.2.data_500more) <- paste("1",colnames(sp.2.data_500more),sep = ".")
colnames(sp.2.data_500more) <- paste("CordBloodCD34",colnames(sp.2.data_500more),sep = "_")
sp.2 <- CreateSeuratObject(counts = Matrix(as.matrix(sp.2.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
sp.2[["percent.mt"]] <- PercentageFeatureSet(sp.2, pattern = "^MT-")
sp.2 <- subset(sp.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
sp.2 <- NormalizeData(sp.2, normalization.method = "LogNormalize", scale.factor = 10000)

sp.normalized <- merge(sp.1, y = sp.2, 
                       add.cell.ids = c("R1", "R2"), project = "Spleen",
                       merge.data = TRUE)
sp.baso <- subset(sp.normalized,subset = IL3RA >0 & FCER1A > 0)
rm(list=setdiff(ls(), "sp.baso"))
##### 
sp.normalized <- FindVariableFeatures(sp.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sp.normalized)
sp.normalized <- ScaleData(sp.normalized, features = all.genes)
sp.normalized <- RunPCA(sp.normalized, features = VariableFeatures(object = sp.normalized))
DimHeatmap(sp.normalized, dims = 1:10, cells = 500, balanced = TRUE)
sp.normalized <- JackStraw(object = sp.normalized, reduction = "pca")
sp.normalized <- ScoreJackStraw(sp.normalized, dims = 1:20)
JackStrawPlot(object = sp.normalized, dims = 1:20)
ElbowPlot(sp.normalized)
sp.normalized <- FindNeighbors(sp.normalized, dims = 1:11)
sp.normalized <- FindClusters(sp.normalized, resolution = 0.5)
sp.normalized <- RunUMAP(sp.normalized, dims = 1:11)
DimPlot(sp.normalized, reduction = "umap")
sp.normalized2 <- JoinLayers(sp.normalized)
sp.normalized2.markers <- FindAllMarkers(sp.normalized2, only.pos = TRUE)

VlnPlot(sp.normalized2, features = c("FCER1A","IL3RA","ITGAM"))

sp.normalized2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
FeaturePlot(sp.normalized, features = c("FCER1A","IL3RA"), blend = TRUE)

DoHeatmap(sp.normalized2, features = top10$gene)

##### Lung #####
lung.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008628_Adult-Lung1_dge.txt.gz",row.names = 1,header = T)
lung.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008629_Adult-Lung2_dge.txt.gz",row.names = 1,header = T)
lung.3.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008630_Adult-Lung3-1_dge.txt.gz",row.names = 1,header = T)
lung.3.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008631_Adult-Lung3-2_dge.txt.gz",row.names = 1,header = T)
lung.3.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008632_Adult-Lung3-3_dge.txt.gz",row.names = 1,header = T)
lung.3.4.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008633_Adult-Lung3-4_dge.txt.gz",row.names = 1,header = T)

lung.1.data_500more = lung.1.data[,colSums(lung.1.data)>=300]
colnames(lung.1.data_500more) <- paste("1",colnames(lung.1.data_500more),sep = ".")
colnames(lung.1.data_500more) <- paste("Lung",colnames(lung.1.data_500more),sep = "_")
lung.1 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.1.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
lung.1[["percent.mt"]] <- PercentageFeatureSet(lung.1, pattern = "^MT-")
lung.1 <- subset(lung.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.1 <- NormalizeData(lung.1, normalization.method = "LogNormalize", scale.factor = 10000)

lung.2.data_500more = lung.2.data[,colSums(lung.2.data)>=300]
colnames(lung.2.data_500more) <- paste("2",colnames(lung.2.data_500more),sep = ".")
colnames(lung.2.data_500more) <- paste("Lung",colnames(lung.2.data_500more),sep = "_")
lung.2 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
lung.2[["percent.mt"]] <- PercentageFeatureSet(lung.2, pattern = "^MT-")
lung.2 <- subset(lung.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.2 <- NormalizeData(lung.2, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.1.data_500more = lung.3.1.data[,colSums(lung.3.1.data)>=300]
colnames(lung.3.1.data_500more) <- paste("3.1",colnames(lung.3.1.data_500more),sep = ".")
colnames(lung.3.1.data_500more) <- paste("Lung",colnames(lung.3.1.data_500more),sep = "_")
lung.3.1 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.1.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.1[["percent.mt"]] <- PercentageFeatureSet(lung.3.1, pattern = "^MT-")
lung.3.1 <- subset(lung.3.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.1 <- NormalizeData(lung.3.1, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.2.data_500more = lung.3.2.data[,colSums(lung.3.2.data)>=300]
colnames(lung.3.2.data_500more) <- paste("3.2",colnames(lung.3.2.data_500more),sep = ".")
colnames(lung.3.2.data_500more) <- paste("Lung",colnames(lung.3.2.data_500more),sep = "_")
lung.3.2 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.2[["percent.mt"]] <- PercentageFeatureSet(lung.3.2, pattern = "^MT-")
lung.3.2 <- subset(lung.3.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.2 <- NormalizeData(lung.3.2, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.3.data_500more = lung.3.3.data[,colSums(lung.3.3.data)>=300]
colnames(lung.3.3.data_500more) <- paste("3.3",colnames(lung.3.3.data_500more),sep = ".")
colnames(lung.3.3.data_500more) <- paste("Lung",colnames(lung.3.3.data_500more),sep = "_")
lung.3.3 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.3.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.3[["percent.mt"]] <- PercentageFeatureSet(lung.3.3, pattern = "^MT-")
lung.3.3 <- subset(lung.3.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.3 <- NormalizeData(lung.3.3, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.4.data_500more = lung.3.4.data[,colSums(lung.3.4.data)>=300]
colnames(lung.3.4.data_500more) <- paste("3.4",colnames(lung.3.4.data_500more),sep = ".")
colnames(lung.3.4.data_500more) <- paste("Lung",colnames(lung.3.4.data_500more),sep = "_")
lung.3.4 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.4.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.4[["percent.mt"]] <- PercentageFeatureSet(lung.3.4, pattern = "^MT-")
lung.3.4 <- subset(lung.3.4, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.4 <- NormalizeData(lung.3.4, normalization.method = "LogNormalize", scale.factor = 10000)

lung.normalized <- merge(lung.1, y = c(lung.2, lung.3.1, lung.3.2, lung.3.3, lung.3.4), 
                         add.cell.ids = c("R1", "R2", "R3", "R4", "R5", "R6"), project = "Lung",
                         merge.data = TRUE)
lung.baso <- subset(lung.normalized,subset = IL3RA >0 & FCER1A > 0)
rm(list=setdiff(ls(), "lung.baso"))
##### 
lung.normalized <- FindVariableFeatures(lung.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(lung.normalized)
lung.normalized <- ScaleData(lung.normalized, features = all.genes)
lung.normalized <- RunPCA(lung.normalized, features = VariableFeatures(object = lung.normalized))
DimHeatmap(lung.normalized, dims = 1:10, cells = 500, balanced = TRUE)
lung.normalized <- JackStraw(object = lung.normalized, reduction = "pca")
lung.normalized <- ScoreJackStraw(lung.normalized, dims = 1:20)
JackStrawPlot(object = lung.normalized, dims = 1:20)
ElbowPlot(lung.normalized)
lung.normalized <- FindNeighbors(lung.normalized, dims = 1:11)
lung.normalized <- FindClusters(lung.normalized, resolution = 0.5)
lung.normalized <- RunUMAP(lung.normalized, dims = 1:11)
DimPlot(lung.normalized, reduction = "umap")
lung.normalized2 <- JoinLayers(lung.normalized)
lung.normalized2.markers <- FindAllMarkers(lung.normalized2, only.pos = TRUE)

VlnPlot(lung.normalized2, features = c("FCER1A","IL3RA","ITGAM"))

lung.normalized2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
FeaturePlot(lung.normalized, features = c("FCER1A","IL3RA"), blend = TRUE)

DoHeatmap(lung.normalized2, features = top10$gene)

baso <- merge(sp.baso, y=c(pbmc.baso,lung.baso,bm.baso),
              project = "Baso",
              merge.data = TRUE)

data.norm <- merge(sp.1, y=c(sp.2,bm.1,bm.2),
                   add.cell.ids = c("SP1","SP2","BM1","BM2"),project="Merge",
                   merge.data = T)

##### Merging all datasets #####
lung.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008628_Adult-Lung1_dge.txt.gz",row.names = 1,header = T)
lung.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008629_Adult-Lung2_dge.txt.gz",row.names = 1,header = T)
lung.3.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008630_Adult-Lung3-1_dge.txt.gz",row.names = 1,header = T)
lung.3.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008631_Adult-Lung3-2_dge.txt.gz",row.names = 1,header = T)
lung.3.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008632_Adult-Lung3-3_dge.txt.gz",row.names = 1,header = T)
lung.3.4.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008633_Adult-Lung3-4_dge.txt.gz",row.names = 1,header = T)
sp.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008649_Adult-Spleen1-1_dge.txt.gz",row.names = 1,header = T)
sp.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008650_Adult-Spleen1-2_dge.txt.gz",row.names = 1,header = T)
bm.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM3943045_Adult-Bone-Marrow1_dge.txt.gz",row.names = 1,header = T)
bm.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM3980128_Adult-Bone-Marrow2_dge.txt.gz",row.names = 1,header = T)
cb34.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008673_Cord-Blood-CD34P1_dge.txt.gz",row.names = 1,header = T)
cb34.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008671_Cord-Blood2-CD34P2-1_dge.txt.gz",row.names = 1,header = T)
cb34.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008672_Cord-Blood2-CD34P2-2_dge.txt.gz",row.names = 1,header = T)
cb.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008668_Cord-Blood1_dge.txt.gz",row.names = 1,header = T)
cb.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008669_Cord-Blood2-1_dge.txt.gz",row.names = 1,header = T)
cb.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008670_Cord-Blood2-2_dge.txt.gz",row.names = 1,header = T)
pbmc.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008638_Adult-Peripheral-Blood1_dge.txt.gz",row.names = 1,header = T)
pbmc.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008639_Adult-Peripheral-Blood2_dge.txt.gz",row.names = 1,header = T)
pbmc.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008640_Adult-Peripheral-Blood3-1_dge.txt.gz",row.names = 1,header = T)
pbmc.3.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008641_Adult-Peripheral-Blood3-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008642_Adult-Peripheral-Blood4-1_dge.txt.gz",row.names = 1,header = T)
pbmc.4.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008643_Adult-Peripheral-Blood4-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008644_Adult-Peripheral-Blood4-3_dge.txt.gz",row.names = 1,header = T)

sp.1.data_500more = sp.1.data[,colSums(sp.1.data)>=300]
colnames(sp.1.data_500more) <- paste("1",colnames(sp.1.data_500more),sep = ".")
colnames(sp.1.data_500more) <- paste("CordBloodCD34",colnames(sp.1.data_500more),sep = "_")
sp.1 <- CreateSeuratObject(counts = Matrix(as.matrix(sp.1.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
sp.1[["percent.mt"]] <- PercentageFeatureSet(sp.1, pattern = "^MT-")
sp.1 <- subset(sp.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
sp.1 <- NormalizeData(sp.1, normalization.method = "LogNormalize", scale.factor = 10000)

sp.2.data_500more = sp.2.data[,colSums(sp.2.data)>=300]
colnames(sp.2.data_500more) <- paste("1",colnames(sp.2.data_500more),sep = ".")
colnames(sp.2.data_500more) <- paste("CordBloodCD34",colnames(sp.2.data_500more),sep = "_")
sp.2 <- CreateSeuratObject(counts = Matrix(as.matrix(sp.2.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
sp.2[["percent.mt"]] <- PercentageFeatureSet(sp.2, pattern = "^MT-")
sp.2 <- subset(sp.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
sp.2 <- NormalizeData(sp.2, normalization.method = "LogNormalize", scale.factor = 10000)

bm.1.data_500more = bm.1.data[,colSums(bm.1.data)>=300]
colnames(bm.1.data_500more) <- paste("1",colnames(bm.1.data_500more),sep = ".")
colnames(bm.1.data_500more) <- paste("CordBloodCD34",colnames(bm.1.data_500more),sep = "_")
bm.1 <- CreateSeuratObject(counts = Matrix(as.matrix(bm.1.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
bm.1[["percent.mt"]] <- PercentageFeatureSet(bm.1, pattern = "^MT-")
bm.1 <- subset(bm.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
bm.1 <- NormalizeData(bm.1, normalization.method = "LogNormalize", scale.factor = 10000)

bm.2.data_500more = bm.2.data[,colSums(bm.2.data)>=300]
colnames(bm.2.data_500more) <- paste("1",colnames(bm.2.data_500more),sep = ".")
colnames(bm.2.data_500more) <- paste("CordBloodCD34",colnames(bm.2.data_500more),sep = "_")
bm.2 <- CreateSeuratObject(counts = Matrix(as.matrix(bm.2.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
bm.2[["percent.mt"]] <- PercentageFeatureSet(bm.2, pattern = "^MT-")
bm.2 <- subset(bm.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
bm.2 <- NormalizeData(bm.2, normalization.method = "LogNormalize", scale.factor = 10000)

cb34.1.data_500more = cb34.1.data[,colSums(cb34.1.data)>=300]
colnames(cb34.1.data_500more) <- paste("1",colnames(cb34.1.data_500more),sep = ".")
colnames(cb34.1.data_500more) <- paste("CordBloodCD34",colnames(cb34.1.data_500more),sep = "_")
cb34.1 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.1.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
cb34.1[["percent.mt"]] <- PercentageFeatureSet(cb34.1, pattern = "^MT-")
cb34.1 <- subset(cb34.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34.1 <- NormalizeData(cb34.1, normalization.method = "LogNormalize", scale.factor = 10000)

cb34.2.data_500more = cb34.2.data[,colSums(cb34.2.data)>=300]
colnames(cb34.2.data_500more) <- paste("1",colnames(cb34.2.data_500more),sep = ".")
colnames(cb34.2.data_500more) <- paste("CordBloodCD34",colnames(cb34.2.data_500more),sep = "_")
cb34.2 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
cb34.2[["percent.mt"]] <- PercentageFeatureSet(cb34.2, pattern = "^MT-")
cb34.2 <- subset(cb34.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34.2 <- NormalizeData(cb34.2, normalization.method = "LogNormalize", scale.factor = 10000)

cb34.3.data_500more = cb34.3.data[,colSums(cb34.3.data)>=300]
colnames(cb34.3.data_500more) <- paste("1",colnames(cb34.3.data_500more),sep = ".")
colnames(cb34.3.data_500more) <- paste("CordBloodCD34",colnames(cb34.3.data_500more),sep = "_")
cb34.3 <- CreateSeuratObject(counts = Matrix(as.matrix(cb34.3.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
cb34.3[["percent.mt"]] <- PercentageFeatureSet(cb34.3, pattern = "^MT-")
cb34.3 <- subset(cb34.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb34.3 <- NormalizeData(cb34.3, normalization.method = "LogNormalize", scale.factor = 10000)

cb.1.data_500more = cb.1.data[,colSums(cb.1.data)>=300]
colnames(cb.1.data_500more) <- paste("1",colnames(cb.1.data_500more),sep = ".")
colnames(cb.1.data_500more) <- paste("CordBloodCD34",colnames(cb.1.data_500more),sep = "_")
cb.1 <- CreateSeuratObject(counts = Matrix(as.matrix(cb.1.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
cb.1[["percent.mt"]] <- PercentageFeatureSet(cb.1, pattern = "^MT-")
cb.1 <- subset(cb.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb.1 <- NormalizeData(cb.1, normalization.method = "LogNormalize", scale.factor = 10000)

cb.2.data_500more = cb.2.data[,colSums(cb.2.data)>=300]
colnames(cb.2.data_500more) <- paste("1",colnames(cb.2.data_500more),sep = ".")
colnames(cb.2.data_500more) <- paste("CordBloodCD34",colnames(cb.2.data_500more),sep = "_")
cb.2 <- CreateSeuratObject(counts = Matrix(as.matrix(cb.2.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
cb.2[["percent.mt"]] <- PercentageFeatureSet(cb.2, pattern = "^MT-")
cb.2 <- subset(cb.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb.2 <- NormalizeData(cb.2, normalization.method = "LogNormalize", scale.factor = 10000)

cb.3.data_500more = cb.3.data[,colSums(cb.3.data)>=300]
colnames(cb.3.data_500more) <- paste("1",colnames(cb.3.data_500more),sep = ".")
colnames(cb.3.data_500more) <- paste("CordBloodCD34",colnames(cb.3.data_500more),sep = "_")
cb.3 <- CreateSeuratObject(counts = Matrix(as.matrix(cb.3.data_500more),sparse=T),
                           min.cells = 3, min.features = 300,names.delim = "\\.")
cb.3[["percent.mt"]] <- PercentageFeatureSet(cb.3, pattern = "^MT-")
cb.3 <- subset(cb.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
cb.3 <- NormalizeData(cb.3, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.1.data_500more = pbmc.1.data[,colSums(pbmc.1.data)>=300]
colnames(pbmc.1.data_500more) <- paste("1",colnames(pbmc.1.data_500more),sep = ".")
colnames(pbmc.1.data_500more) <- paste("CordBloodCD34",colnames(pbmc.1.data_500more),sep = "_")
pbmc.1 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.1.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.1[["percent.mt"]] <- PercentageFeatureSet(pbmc.1, pattern = "^MT-")
pbmc.1 <- subset(pbmc.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.1 <- NormalizeData(pbmc.1, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.2.data_500more = pbmc.2.data[,colSums(pbmc.2.data)>=300]
colnames(pbmc.2.data_500more) <- paste("1",colnames(pbmc.2.data_500more),sep = ".")
colnames(pbmc.2.data_500more) <- paste("CordBloodCD34",colnames(pbmc.2.data_500more),sep = "_")
pbmc.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.2, pattern = "^MT-")
pbmc.2 <- subset(pbmc.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.2 <- NormalizeData(pbmc.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.3.data_500more = pbmc.3.data[,colSums(pbmc.3.data)>=300]
colnames(pbmc.3.data_500more) <- paste("1",colnames(pbmc.3.data_500more),sep = ".")
colnames(pbmc.3.data_500more) <- paste("CordBloodCD34",colnames(pbmc.3.data_500more),sep = "_")
pbmc.3 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.3.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.3[["percent.mt"]] <- PercentageFeatureSet(pbmc.3, pattern = "^MT-")
pbmc.3 <- subset(pbmc.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.3 <- NormalizeData(pbmc.3, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.3.2.data_500more = pbmc.3.2.data[,colSums(pbmc.3.2.data)>=300]
colnames(pbmc.3.2.data_500more) <- paste("1",colnames(pbmc.3.2.data_500more),sep = ".")
colnames(pbmc.3.2.data_500more) <- paste("CordBloodCD34",colnames(pbmc.3.2.data_500more),sep = "_")
pbmc.3.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.3.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.3.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.3.2, pattern = "^MT-")
pbmc.3.2 <- subset(pbmc.3.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.3.2 <- NormalizeData(pbmc.3.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.data_500more = pbmc.4.data[,colSums(pbmc.4.data)>=300]
colnames(pbmc.4.data_500more) <- paste("1",colnames(pbmc.4.data_500more),sep = ".")
colnames(pbmc.4.data_500more) <- paste("CordBloodCD34",colnames(pbmc.4.data_500more),sep = "_")
pbmc.4 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4[["percent.mt"]] <- PercentageFeatureSet(pbmc.4, pattern = "^MT-")
pbmc.4 <- subset(pbmc.4, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4 <- NormalizeData(pbmc.4, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.2.data_500more = pbmc.4.2.data[,colSums(pbmc.4.2.data)>=300]
colnames(pbmc.4.2.data_500more) <- paste("1",colnames(pbmc.4.2.data_500more),sep = ".")
colnames(pbmc.4.2.data_500more) <- paste("CordBloodCD34",colnames(pbmc.4.2.data_500more),sep = "_")
pbmc.4.2 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4.2[["percent.mt"]] <- PercentageFeatureSet(pbmc.4.2, pattern = "^MT-")
pbmc.4.2 <- subset(pbmc.4.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4.2 <- NormalizeData(pbmc.4.2, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.4.3.data_500more = pbmc.4.3.data[,colSums(pbmc.4.3.data)>=300]
colnames(pbmc.4.3.data_500more) <- paste("1",colnames(pbmc.4.3.data_500more),sep = ".")
colnames(pbmc.4.3.data_500more) <- paste("CordBloodCD34",colnames(pbmc.4.3.data_500more),sep = "_")
pbmc.4.3 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.4.3.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
pbmc.4.3[["percent.mt"]] <- PercentageFeatureSet(pbmc.4.3, pattern = "^MT-")
pbmc.4.3 <- subset(pbmc.4.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
pbmc.4.3 <- NormalizeData(pbmc.4.3, normalization.method = "LogNormalize", scale.factor = 10000)

lung.1.data_500more = lung.1.data[,colSums(lung.1.data)>=300]
colnames(lung.1.data_500more) <- paste("1",colnames(lung.1.data_500more),sep = ".")
colnames(lung.1.data_500more) <- paste("CordBloodCD34",colnames(lung.1.data_500more),sep = "_")
lung.1 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.1.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
lung.1[["percent.mt"]] <- PercentageFeatureSet(lung.1, pattern = "^MT-")
lung.1 <- subset(lung.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.1 <- NormalizeData(lung.1, normalization.method = "LogNormalize", scale.factor = 10000)

lung.2.data_500more = lung.2.data[,colSums(lung.2.data)>=300]
colnames(lung.2.data_500more) <- paste("1",colnames(lung.2.data_500more),sep = ".")
colnames(lung.2.data_500more) <- paste("CordBloodCD34",colnames(lung.2.data_500more),sep = "_")
lung.2 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.2.data_500more),sparse=T),
                             min.cells = 3, min.features = 300,names.delim = "\\.")
lung.2[["percent.mt"]] <- PercentageFeatureSet(lung.2, pattern = "^MT-")
lung.2 <- subset(lung.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.2 <- NormalizeData(lung.2, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.1.data_500more = lung.3.1.data[,colSums(lung.3.1.data)>=300]
colnames(lung.3.1.data_500more) <- paste("1",colnames(lung.3.1.data_500more),sep = ".")
colnames(lung.3.1.data_500more) <- paste("CordBloodCD34",colnames(lung.3.1.data_500more),sep = "_")
lung.3.1 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.1.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.1[["percent.mt"]] <- PercentageFeatureSet(lung.3.1, pattern = "^MT-")
lung.3.1 <- subset(lung.3.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.1 <- NormalizeData(lung.3.1, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.2.data_500more = lung.3.2.data[,colSums(lung.3.2.data)>=300]
colnames(lung.3.2.data_500more) <- paste("1",colnames(lung.3.2.data_500more),sep = ".")
colnames(lung.3.2.data_500more) <- paste("CordBloodCD34",colnames(lung.3.2.data_500more),sep = "_")
lung.3.2 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.2.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.2[["percent.mt"]] <- PercentageFeatureSet(lung.3.2, pattern = "^MT-")
lung.3.2 <- subset(lung.3.2, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.2 <- NormalizeData(lung.3.2, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.3.data_500more = lung.3.3.data[,colSums(lung.3.3.data)>=300]
colnames(lung.3.3.data_500more) <- paste("1",colnames(lung.3.3.data_500more),sep = ".")
colnames(lung.3.3.data_500more) <- paste("CordBloodCD34",colnames(lung.3.3.data_500more),sep = "_")
lung.3.3 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.3.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.3[["percent.mt"]] <- PercentageFeatureSet(lung.3.3, pattern = "^MT-")
lung.3.3 <- subset(lung.3.3, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.3 <- NormalizeData(lung.3.3, normalization.method = "LogNormalize", scale.factor = 10000)

lung.3.4.data_500more = lung.3.4.data[,colSums(lung.3.4.data)>=300]
colnames(lung.3.4.data_500more) <- paste("1",colnames(lung.3.4.data_500more),sep = ".")
colnames(lung.3.4.data_500more) <- paste("CordBloodCD34",colnames(lung.3.4.data_500more),sep = "_")
lung.3.4 <- CreateSeuratObject(counts = Matrix(as.matrix(lung.3.4.data_500more),sparse=T),
                               min.cells = 3, min.features = 300,names.delim = "\\.")
lung.3.4[["percent.mt"]] <- PercentageFeatureSet(lung.3.4, pattern = "^MT-")
lung.3.4 <- subset(lung.3.4, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
lung.3.4 <- NormalizeData(lung.3.4, normalization.method = "LogNormalize", scale.factor = 10000)



all.data <- merge(sp.1, y=c(sp.2,bm.1,bm.2,cb34.1,cb34.2,cb34.3,cb.1,cb.2,cb.3,
                            lung.1,lung.2,lung.3.1,lung.3.2,lung.3.3,lung.3.4,
                            pbmc.1,pbmc.2,pbmc.3,pbmc.3.2,pbmc.4,pbmc.4.2,pbmc.4.3),
                  add.cell.ids=c("Spleen1","Spleen2","BoneMarrow1","BoneMarrow2",
                                 "CB34_1","CB34_2","CB34_3","CB1","CB2","CB3",
                                 "Lung1","Lung2","Lung3_1","Lung3_2","Lung3_3","Lung3_4",
                                 "PBMC1","PBMC2","PBMC3","PBMC3_2","PBMC4","PBMC4_2","PBMC4_3"),
                  project="Baso",
                  merge.data = T)
all.data <- JoinLayers(all.data)
all.data.baso <- subset(all.data, subset = IL3RA >0 & FCER1A > 0)
all.data.baso <- FindVariableFeatures(all.data.baso, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all.data)
all.data.baso <- ScaleData(all.data.baso, features = all.genes)
all.data.baso <- RunPCA(all.data.baso, features = VariableFeatures(object = all.data.baso))
DimHeatmap(all.data.baso, dims = 1:10, cells = 500, balanced = TRUE)
all.data <- JackStraw(object = all.data, reduction = "pca")
all.data <- ScoreJackStraw(all.data, dims = 1:20)
JackStrawPlot(object = all.data, dims = 1:20)
ElbowPlot(all.data)
all.data.baso <- FindNeighbors(all.data.baso, dims = 1:11)
all.data.baso <- FindClusters(all.data.baso, resolution = 0.5)
all.data.baso <- RunUMAP(all.data.baso, dims = 1:11)
DimPlot(all.data.baso, reduction = "umap")
all.data2 <- JoinLayers(all.data)
all.data2.markers <- FindAllMarkers(all.data2, only.pos = TRUE)