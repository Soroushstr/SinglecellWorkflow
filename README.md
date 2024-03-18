# SinglecellWorkflow
This repository is a notebook for me for establishment of a proper single-cell RNASeq data analysis.

# Installing and loading libraries
```R
install.packages("Seurat")
devtools::install_github('immunogenomics/presto')
devtools::install_github("ImmuneDynamics/Spectre")
library(Seurat)
library(Matrix)
library(dplyr)
```

# Loading the data 
All peripheral blood datasets were preprocessed and normalized and then merged for further analysis and establishment of a reference Basophil dataset. 
```R
pbmc.1.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008638_Adult-Peripheral-Blood1_dge.txt.gz",row.names = 1,header = T)
```
# Preprocessing
## QC
Cells with at least 300 Unique Molecular Identifiers (UMIs) were selected from each peripheral blood dataset.
```R
pbmc.1.data_500more = pbmc.1.data[,colSums(pbmc.1.data)>=300]
colnames(pbmc.1.data_500more) <- paste("1",colnames(pbmc.1.data_500more),sep = ".")
colnames(pbmc.1.data_500more) <- paste("PeriBlood",colnames(pbmc.1.data_500more),sep = "_")
```

Cells with at least 300 and lower than 2500 transcripts, and transcripts present in at least 3 cells were selected as a subset of each dataset. 
```R
pbmc.1 <- CreateSeuratObject(counts = Matrix(as.matrix(pbmc.1.data_500more),sparse=T),
                            min.cells = 3, min.features = 300,names.delim = "\\.")
```
## Removing cells with high mitochondrial content
Creating a column corresponding to mitochondrial gene content to metadata
Selecting cells with <20% mitochondrial content
```R
pbmc.1[["percent.mt"]] <- PercentageFeatureSet(pbmc.1, pattern = "^MT-")
pbmc.1 <- subset(pbmc.1, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 20)
```

## Normalization
Log-normalizing counts
```R
pbmc.1 <- NormalizeData(pbmc.1, normalization.method = "LogNormalize", scale.factor = 10000)
```

## Merging
If you are using only one dataset, you can only do following line of code and continue your analysis
```R
pbmc.normalized <- pbmc.1
```
Once the preprocessing and normalization was done for all datasets, they get merged into one dataset
```R
pbmc.normalized <- merge(pbmc.1, y = c(pbmc.2, pbmc.3, pbmc.3.2, pbmc.4, pbmc.4.2, pbmc.4.3), 
                         add.cell.ids = c("R1", "R2", "R3", "R3.2", "R4", "R4.2", "R4.3"), project = "PBMC12K",
                         merge.data = TRUE)
all.genes <- rownames(pbmc.normalized)
```

## Finding Variable Genes
To clarify difference between cell types, top 2000 variable genes, based on the feature variances are selected for downstream analyses
```R
pbmc.normalized <- FindVariableFeatures(pbmc.normalized, selection.method = "vst", nfeatures = 2000)
```

## Scale and center data
```R
pbmc.normalized <- ScaleData(pbmc.normalized, features = all.genes)
```

# Analysis and Processing
## PCA
Doing Principle Component Analysis (PCA) and plotting heatmap of top 10 PCs
```R
pbmc.normalized <- RunPCA(pbmc.normalized, features = VariableFeatures(object = pbmc.normalized))
DimHeatmap(pbmc.normalized, dims = 1:10, cells = 500, balanced = TRUE)
```

## Selecting PCs
To select significant Principle Components (PCs), Jackstraw Plot and Elbow Plot methods are utilized.
### Jackstraw method
```R
pbmc.normalized <- JackStraw(object = pbmc.normalized, reduction = "pca")
pbmc.normalized <- ScoreJackStraw(pbmc.normalized, dims = 1:20)
JackStrawPlot(object = pbmc.normalized, dims = 1:20)
```
Jackstraw plot of the Principle Components (PCs) of the peripheral samples:

!(JackstrawPeripheralBloodPCs)[https://github.com/Soroushstr/SinglecellWorkflow/blob/ecf915e63fc0135155dfdb4e51e4b98e8b8e8a15/PBMC%20Jackstraw.png]
### Elbow plot method
```R
ElbowPlot(pbmc.normalized)
```
## Clustering
Finding neighbors and clustering cells
```R
pbmc.normalized <- FindNeighbors(pbmc.normalized, dims = 1:11)
pbmc.normalized <- FindClusters(pbmc.normalized, resolution = 0.5)
```
## UMAP
Running UMAP and Plotting cells in the UMAP representation
```R
pbmc.normalized <- RunUMAP(pbmc.normalized, dims = 1:11)
DimPlot(pbmc.normalized, reduction = "umap")
```
![UMAPPBMC](https://github.com/Soroushstr/SinglecellWorkflow/blob/4135c99ad056b70f236e29c78377d0380bf0ef98/UMAP%20PBMC.png)

# Processing single-cell data
## Lymph nodes data (Kim et al)

### Installing and importing required packages
```Python
! pip install scanpy
! pip install bbknn
! pip install gprofiler
! pip install sfaira
```
```Python
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import bbknn
import anndata
from gprofiler import gprofiler

import warnings

from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects import pandas2ri

%load_ext rpy2.ipython

warnings.filterwarnings("ignore", category=RRuntimeWarning)
pandas2ri.activate()

plt.rcParams['figure.figsize']=(8,8)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
sc.settings.verbosity = 3
```
