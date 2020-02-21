#date created: Monday, 16/12/2019
#author: Maren Buettner (ICB)
#project: Early life stress in collaboration with Pablo Lopez (Chenlab)
#purpose: We compute size factors for all samples based on a joint louvain clustering. These clusters were created in
# 01_preprocessing.ipynb
# 


library(rhdf5)
library(Matrix)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scran)

#set paths and dates
figure_path <- '~/Documents/Collaborations/Lopez_Chenlab/figures/'
file_path <- '~/Documents/Collaborations/Lopez_Chenlab/table/'
data_path <- "~/Documents/Collaborations/Lopez_Chenlab/data/"
today <- format(Sys.Date(), '%y%m%d')

###################
# load data files #
###################

#final annotation data 
f_path <- paste0(data_path, 'tmp.h5ad')
adata <- h5read(f_path,'/',compoundAsDataFrame=FALSE)

#get genes/cell ID
barcodes <- unlist(adata$obs$index)
genes <- unlist(adata$var$index)

data_counts <- adata$X$data
index_counts <- adata$X$indices
ptr_counts <- adata$X$indptr
sparse_mat <- sparseMatrix(p = as.numeric(ptr_counts), 
                           x = as.numeric(data_counts), 
                           i = as.numeric(index_counts)+1)

louvain <- adata$obs$louvain

#compute n_genes, n_counts
n_genes <- as.numeric(adata$obs$n_genes)
n_counts <- as.numeric(adata$obs$n_counts)

#compute size factors with scran
size_factors = computeSumFactors(as.matrix(sparse_mat), clusters=louvain, min.mean=0.1)

#get attributes
cellData <- data.frame(
  sample=factor(unlist(adata$uns$sample_name_categories)[unlist(adata$obs$sample_name)+1]),
  condition = factor(unlist(adata$uns$condition_categories)[unlist(adata$obs$condition)+1]),
  n_genes = n_genes,
  n_counts = n_counts,
  #size_factors = size_factors,
  #cell_type = factor(unlist(adata$uns$cell_type_categories)[unlist(adata$obs$cell_type+1)]),
  louvain = factor(unlist(adata$uns$louvain_categories)[unlist(adata$obs$louvain+1)])
) 

ggplot(cellData, aes(louvain, size_factors, fill=sample)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#save size factors to file
write.csv(x = size_factors, file = paste0(data_path, 'size_factors_all_samples.csv'), quote = FALSE, row.names = barcodes)
