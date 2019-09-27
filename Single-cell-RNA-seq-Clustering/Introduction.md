# Single-cell-RNA-seq-Clustering (R/matlab/python)

Order: 1.Corr → 2.MPSSC → 3.NMF → 4.PhenoGraph → 5.RAFSIL → 6.SC3 → 7.Seurat → 8.SIMLR → 9.SNN-cliq

## 1. Corr (matlab)

Ref: Jiang H, Sohn L L, Huang H, et al. Single cell clustering based on cell-pair differentiability correlation and variance analysis[J]. Bioinformatics, 2018, 34(21): 3684-3694.

```linux

%% Example using Biase data set:
clear all

%% load data sets ('in_X' and 'Labels') (./Corr/Demo.mat)
load('Demo.mat')

%% Note: one can use any data set that consists of in_X and true_labs, where in_X is an n by p gene expression matrix and true_labs is the ground truth labels. Here n and p are number of cells and genes, respectively.

%% Automatically determine the optimal number of cluster(slowly)
optnumber=Number_Corr(in_X);
cluster = optnumber;

%% or setup by user
cluster = *;

%% Clustering
clusters = Corr_Clustering(in_X,cluster,Labels);

%% evaluation(./Corr/Cal_NMI.m; ./Corr/RandIndex.m)
NMI = Cal_NMI(true_labs,clusters);
ARI = RandIndex(true_labs,clusters);
```


## 2. MPSSC (matlab)

Ref: Park S, Zhao H. Spectral clustering based on learning similarity matrix[J]. Bioinformatics, 2018, 34(12): 2069-2076.

```linux

%% Example using Biase data set:
clear all
addpath(genpath(pwd))

%% load data sets ('in_X' and 'true_labs') (./MPSSC/Biase.mat)
load('Biase.mat')

%% Note: one can use any data set that consists of in_X and true_labs, where in_X is an n by p gene expression matrix and true_labs is the ground truth labels. Here n and p are number of cells and genes, respectively.

%% Penalty parameters:
rho=0.2; lam=0.0001; lam2=lam; eta=1; c=0.1;

%% Run MPSSC and obtain the target matrix P
[P] = clus_fin_update(rho, lam, lam2, eta, c, in_X, true_labs); 

%% Obtain clustering labels "clus_labs" and compute NMI measure:
[NMI, ~,clus_labs,~]=calc2_nmis(max(true_labs), double(P),true_labs);   

%% Compute performance measures
%% Compute Purity
Purity=purity(max(true_labs), clus_labs, true_labs)

%% Compute ARI
ARI=RandIndex(clus_labs,true_labs)

%%% Final output: Performances (three measures) of MPSSC
[NMI, ARI]
```


## 3. NMF (python & R)

Ref: Shao C, Höfer T. Robust classification of single-cell transcriptome data by nonnegative matrix factorization[J]. Bioinformatics, 2017, 33(2): 235-242. 

```linux

## Example using Biase data set:

## Code
unzip file nimfa-master.zip(./NMF/nimfa-master.zip)

## load data
Copy Demo data (./NMF/Data_Biase_input_.csv; ./NMF/Data_Biase_labels.csv) to the nimfa-master\sc-RNAseq_scripts folder

## Note: Multiple data can be copied, and the method runs accordingly.

## Run Command:
python NMF_scRNAseq.py --min_rank=7 --max_rank=7 --max_depth=1 --RHeatmap=".\sc-RNAseq_scripts\nmfHeatmap.R"

## Note: The first few parameters are automatically calculated and will be automatically extracted from the data and labels. The final folder path should be set to the location of the nmfHeatmap.R file.

## Final output: *(data name) depth_*_r*.label_input_.csv format:
		cells	cls
sample1	sample1	3
sample2	sample2	5
sample3	sample3	3
## The third column is the classification result.
```


## 4. PhenoGraph (Python)

Ref: Levine J H, Simonds E F, Bendall S C, et al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis[J]. Cell, 2015, 162(1): 184-197.

```linux
## install PhenoGraph, simply run the setup script:
python3 setup.py install

## Or use:
pip3 install git+https://github.com/jacoblevine/phenograph.git

## Expected use is within a script or interactive kernel running Python `3.x`. Data are expected to be passed as a `numpy.ndarray`. When applicable, the code uses CPU multicore parallelism via `multiprocessing`. 
## To run basic clustering (data: ./PhenoGraph/Biase.csv):
import numpy as np
import phenograph

if __name__ == "__main__":
    data = np.loadtxt(open("Biase.csv", "r"), delimiter=",")
    data = data.T
    communities, graph, Q = phenograph.cluster(data)
	print(communities)

## For a dataset of *N* rows, `communities` will be a length *N* vector of integers specifying a community assignment for each row in the data. Any rows assigned `-1` were identified as *outliers* and should not be considered as a member of any community.
## `graph` is a *N* x *N* `scipy.sparse` matrix representing the weighted graph used for community detection. 
## `Q` is the modularity score for `communities` as applied to `graph`.
```


## 5. RAFSIL (R)

Ref: Pouyan M B, Kostka D. Random forest based similarity learning for single cell RNA sequencing data[J]. Bioinformatics, 2018, 34(13): i79-i88.

```linux
## data: an (n x m) data matrix of gene expression measurements of single cells
## NumC: number of clusters to be extracted over data
## nrep: integer for RAFSIL1 the number of times similarities are learned
## method: character "RAFISIL1" or "RAFSIL2"
## gene_filter: logical should gene filtering be performed?
## load data (./RAFSIL/Biase.Rdata)
load("Biase.Rdata")

## Clustering
library(RAFSIL)
res1 = RAFSIL(data, NumC = NULL, method = "RAFSIL1", gene_filter = TRUE)
res2 = RAFSIL(data, NumC = NULL, method = "RAFSIL2", gene_filter = TRUE)

## output results
## note: if NumC==NULL there no clustering results
## dissimilarity matrix
res1$D

## Clustering Results
res1$lab

## evaluation NMI, RI, ARI(./evalcluster.R)
## evalcluster(truelabel,predlabel) 
evalcluster(label,res1$lab)
```


## 6. SC3 (R)

Ref: Kiselev V Y, Kirschner K, Schaub M T, et al. SC3: consensus clustering of single-cell RNA-seq data[J]. Nature methods, 2017, 14(5): 483.

```linux
## install SC3
source("https://bioconductor.org/biocLite.R")
biocLite("SC3")

## load data (./SC3/Demo.Rdata)
load(Demo.Rdata)

## create a SingleCellExperiment object
library(SC3)
library(SingleCellExperiment)
library(scater)

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(yan),
        logcounts = log2(as.matrix(yan) + 1)
    ), 
    colData = ann
)

## define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)

## remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

## Clustering
sce <- sc3(sce, ks = *,gene_filter = FALSE,biology = FALSE)

## note that ks is the number of clusters

## cluster results (* is the number of cluster which seted by user)
colData(sce)$sc3_*_clusters

## evaluation NMI, RI, ARI(./evalcluster.R)
## evalcluster(truelabel,predlabel) 
evalcluster(label,colData(sce)$sc3_*_clusters)
```


## 7. Seurat (R)

Ref: Satija R, Farrell J A, Gennert D, et al. Spatial reconstruction of single-cell gene expression data[J]. Nature biotechnology, 2015, 33(5): 495.

```linux
library(Seurat)
library(dplyr)

## load data (./Seurat/Koh.Rdata)
## Unnormalized data such as raw counts or TPMs
dataname = "Koh.Rdata"
load(dataname)

## create Seurat object
pbmc_small <- CreateSeuratObject(data)

## Normalize the count data present in a given assay.
pbmc_small <- NormalizeData(object = pbmc_small)

## Identifies features that are outliers on a 'mean variability plot'.
pbmc_small <- FindVariableFeatures(object = pbmc_small)

## Scales and centers features in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each feautre, and the resulting residuals are then scaled and centered.
pbmc_small <- ScaleData(
  object = pbmc_small
)

## Run a PCA dimensionality reduction. For details about stored PCA calculation parameters, see PrintPCAParams.
pbmc_small <- RunPCA(
  object = pbmc_small,
  pc.genes = pbmc_small@var.genes)

## Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component.
pbmc_small <- JackStraw(
  object = pbmc_small)

## Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
pbmc_small <- FindNeighbors(pbmc_small)

##Clustering
res = FindClusters(object = pbmc_small)
res$seurat_clusters

## evaluation NMI, RI, ARI(./evalcluster.R)
## evalcluster(truelabel,predlabel) 
evalcluster(label,res$seurat_clusters)
```


## 8. SIMLR (R)

Ref: Wang B, Zhu J, Pierson E, et al. Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning[J]. Nature methods, 2017, 14(4): 414.

```linux
## install
source("https://bioconductor.org/biocLite.R")
biocLite("SIMLR")

## load data (./SIMLR/Biase.Rdata)
load("Biase.Rdata")

## Run
library(SIMLR)
library(igraph)

set.seed(11111)
example = SIMLR(X = data, c = length(table(label)),normalize = TRUE)

## cluster results
example$y$cluster

## similarity computed by SIMLR
example$S

# plot 
plot(example$ydata,
     col = c(topo.colors(k))[label],
     xlab = "SIMLR component 1",
     ylab = "SIMLR component 2",
     pch = 20,
     main="SIMILR 2D visualization for data")

## evaluation NMI, RI, ARI(./evalcluster.R)
## evalcluster(truelabel,predlabel) 
evalcluster(label,example$y$cluster)
```



## 9. SNN-cliq (R & python)

Ref: Xu C, Su Z. Identification of cell types from single-cell transcriptomes using a novel clustering method[J]. Bioinformatics, 2015, 31(12): 1974-1980.

```linux
## load data (./SNN-cliq/Biase.Rdata)
data = load("Biase.Rdata")

## construct Edge file(KNN graph)
# ./SNN-cliq/SNN.R
source('SNN.R')
# edge_file is the output Edge file name, such as "./BiaseEdge"
SNN(data, edge_file, k=3, distance='euclidean')

## Clusterig by Cliq method (./SNN-cliq/Cliq.py)
python Cliq.py -i edge_file -o <out_file> -r 0.7 -m 0.5
# note that 
# usage: SNNgraph_clustering.py -i <edge_file> -o <out_file> [options]
#  -i,--input    input file path and name
#  -o,--output   output file path and name
#  optional arguments:
#  -r,--r-quasi-cliq     quasi-clique parameter. A number in the range of (0  1]. Default is 0.7.
#  -m,--merging  merging parameter. A number in the range of (0  1]. Default is 0.5.
#  -n,--number   number of objects. Default is the maximum index in the input file.
#  -h,--help     print help message.
```
