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
pbmc.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008639_Adult-Peripheral-Blood2_dge.txt.gz",row.names = 1,header = T)
pbmc.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008640_Adult-Peripheral-Blood3-1_dge.txt.gz",row.names = 1,header = T)
pbmc.3.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008641_Adult-Peripheral-Blood3-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008642_Adult-Peripheral-Blood4-1_dge.txt.gz",row.names = 1,header = T)
pbmc.4.2.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008643_Adult-Peripheral-Blood4-2_dge.txt.gz",row.names = 1,header = T)
pbmc.4.3.data <- read.table(file = "C:/Users/sh/Downloads/Old Laptop/EPFL/Baso single cell/GSE134335/GSM4008644_Adult-Peripheral-Blood4-3_dge.txt.gz",row.names = 1,header = T)
```
# QC
```R
# Selecting cells with at least 500 transcripts
adata.data_500more = adata.data[,colSums(adata.data)>=500]
colnames(adata.data_500more) <- paste("2",colnames(adata.data_500more),sep = ".")
colnames(adata.data_500more) <- paste("CordBlood",colnames(adata.data_500more),sep = "_")
adata <- CreateSeuratObject(counts = Matrix(as.matrix(adata.data_500more),sparse=T),
                            min.cells = 3, min.features = 300,names.delim = "\\.")
# Adding Mitochondrial percentage metadata
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
# QC plots
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```
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
