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
