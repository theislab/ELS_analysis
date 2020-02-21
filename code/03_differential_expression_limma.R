library(rhdf5)
library(Matrix)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#set paths and dates
figure_path <- '~/Documents/Collaborations/Lopez_Chenlab/figures/'
file_path <- '~/Documents/Collaborations/Lopez_Chenlab/table/'
data_path <- "~/Documents/Collaborations/Lopez_Chenlab/data/"
today <- format(Sys.Date(), '%y%m%d')

###################
# load data files #
###################

#final annotation data 
f_path <- paste0(data_path, 'cellxgene.h5ad')
adata <- h5read(f_path,'/',compoundAsDataFrame=FALSE)

#read size factors (computed on full sample)
size_factors <- read.csv(paste0(data_path, 'size_factors_all_samples.csv'))$x

#get genes/cell ID
barcodes <- unlist(adata$obs$index)
genes <- unlist(adata$var$index)

data_counts <- adata$layers$counts$data
index_counts <- adata$layers$counts$indices
ptr_counts <- adata$layers$counts$indptr
sparse_mat <- sparseMatrix(p = as.numeric(ptr_counts), 
                           x=as.numeric(data_counts), 
                           i = as.numeric(index_counts)+1)

#size_factors <- adata$obs$size_factors
louvain <- adata$obs$louvain

#compute n_genes, n_counts
n_genes <- as.numeric(adata$obs$n_genes)
n_counts <- as.numeric(adata$obs$n_counts)

#set levels
#sample_levels <- c( "mRFP","mTmG", "mGFP")
#cell_type_levels <- c("alpha", "beta", "gamma", "delta" , "epsilon",
#                      "alpha-delta", "beta-delta", "alpha-beta",          
#                      "alpha, proliferating", "beta, proliferating",
#                      "gamma, proliferating", "delta, proliferating")

#get attributes
cellData <- data.frame(
  sample=factor(unlist(adata$uns$sample_name_categories)[unlist(adata$obs$sample_name)+1]),
  baseline=factor(unlist(adata$uns$baseline_categories)[unlist(adata$obs$baseline)+1]),
  adult_stress=factor(unlist(adata$uns$adult_stress_categories)[unlist(adata$obs$adult_stress)+1]),
  condition = factor(unlist(adata$uns$condition_categories)[unlist(adata$obs$condition)+1]),
  n_genes = n_genes,
  n_counts = n_counts,
  size_factors = size_factors,
  cell_type = factor(unlist(adata$uns$cell_type_categories)[unlist(adata$obs$cell_type+1)]),
  louvain = factor(unlist(adata$uns$louvain_categories)[unlist(adata$obs$louvain+1)])
) 




#create boxplots of genes and library size per cell type
ggplot(cellData, aes(cell_type, n_genes, fill=sample)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_genes, fill=baseline)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_genes, fill=adult_stress)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_genes, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, size_factors, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(louvain, n_counts, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_counts, fill=sample)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_counts, fill=baseline)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_counts, fill=adult_stress)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData, aes(cell_type, n_counts, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#############################
#pairwise test of all cells #
#############################
logFC_thresh <- 0.05
p_adj_thresh <- 1e-5

#prepare pairwise tests
celltypes_OI <- levels(cellData$cell_type)
#select by baseline

for (baseline_idx in levels(cellData$baseline)){
  print(paste0('Differential expression test on ', baseline_idx))
  tmp_idx <- cellData$baseline %in% baseline_idx
  cellData_tmp <- cellData[tmp_idx,]
  #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
  sparse_mat_tmp <- sparse_mat[,tmp_idx]
  
  for (i in celltypes_OI){
    cell_types <-i
    print(paste0('Test in ', cell_types, ' cells.'))
    
    cell_idx <- cellData_tmp$cell_type %in% cell_types
    data_test <- as.matrix(sparse_mat_tmp[,cell_idx])
    #scale data to log-scran
    data_test <- log(data_test %*% diag(1/cellData_tmp$size_factors[cell_idx]) + 1)
    rownames(data_test) <- genes
    
    keep_genes <- rowMeans(data_test)>0.05
    print(paste0('Keep ', sum(keep_genes), ' genes.'))
    data_test <- data_test[keep_genes,]  
    
    condition_test <- cellData_tmp$adult_stress[cellData_tmp$cell_type %in% cell_types]
    condition_test <- droplevels(condition_test)
    
    rep_test <- cellData_tmp$sample[cellData_tmp$cell_type %in% cell_types]
    rep_test <- droplevels(rep_test)
    #ct_test <- cellData$progen_test_type[cellData_tmp$progen_test_type %in% cell_types]
    #ct_test <- droplevels(ct_test)
    #cell_types_test <- levels(ct_test)
    #levels(ct_test) <- c('ct1', 'ct2')
    #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
    #?limma
    #prepare model matrix
    
    dge <- edgeR::DGEList(data_test,sample=rep_test, group = condition_test)
    design <- model.matrix(~condition_test + rep_test)
    colnames(design) <- gsub('condition_test','', colnames(design))
    #contrasts <- makeContrasts(ct1-ct2, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef= 2,
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>logFC_thresh & rest$adj.P.Val<p_adj_thresh,]
    title <- paste0(baseline_idx, '_', colnames(fit$coefficients)[2],'_vs_', 
                    names(table(condition_test))[!names(table(condition_test))==colnames(fit$coefficients)[2]], 
                    '_', cell_types
    )
    write.csv(x = subset_rest, 
              file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
    
    #create volcano plot
    data_res <- subset_rest %>% select(logFC, adj.P.Val)
    data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
    
    g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
      labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
      theme_classic()
    ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
           plot = g0, width = 6, height=6)
  }
}

#select by adult_stress

for (stress_idx in levels(cellData$adult_stress)){
  print(paste0('Differential expression test on ', stress_idx))
  tmp_idx <- cellData$adult_stress %in% stress_idx
  cellData_tmp <- cellData[tmp_idx,]
  #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
  sparse_mat_tmp <- sparse_mat[,tmp_idx]
  
  for (i in celltypes_OI){
    cell_types <-i
    print(paste0('Test in ', cell_types, ' cells.'))
    
    cell_idx <- cellData_tmp$cell_type %in% cell_types
    data_test <- as.matrix(sparse_mat_tmp[,cell_idx])
    #scale data to log-scran
    data_test <- log(data_test %*% diag(1/cellData_tmp$size_factors[cell_idx]) + 1)
    rownames(data_test) <- genes
    
    keep_genes <- rowMeans(data_test)>0.05
    print(paste0('Keep ', sum(keep_genes), ' genes.'))
    data_test <- data_test[keep_genes,]  
    
    condition_test <- cellData_tmp$baseline[cellData_tmp$cell_type %in% cell_types]
    condition_test <- droplevels(condition_test)
    
    rep_test <- cellData_tmp$sample[cellData_tmp$cell_type %in% cell_types]
    rep_test <- droplevels(rep_test)
    #ct_test <- cellData$progen_test_type[cellData_tmp$progen_test_type %in% cell_types]
    #ct_test <- droplevels(ct_test)
    #cell_types_test <- levels(ct_test)
    #levels(ct_test) <- c('ct1', 'ct2')
    #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
    #?limma
    #prepare model matrix
    
    dge <- edgeR::DGEList(data_test,sample=rep_test, group = condition_test)
    design <- model.matrix(~condition_test + rep_test)
    colnames(design) <- gsub('condition_test','', colnames(design))
    #contrasts <- makeContrasts(ct1-ct2, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef= 2,
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>logFC_thresh & rest$adj.P.Val<p_adj_thresh,]
    title <- paste0(stress_idx, '_', colnames(fit$coefficients)[2],'_vs_', 
                    names(table(condition_test))[!names(table(condition_test))==colnames(fit$coefficients)[2]], 
                    '_', cell_types
    )
    write.csv(x = subset_rest, 
              file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
    
    #create volcano plot
    data_res <- subset_rest %>% select(logFC, adj.P.Val)
    data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
    
    g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
      labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
      theme_classic()
    ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
           plot = g0, width = 6, height=6)
  }
}

##############################
# test glutamatergic neurons #
##############################
logFC_thresh <- 0.05
p_adj_thresh <- 1e-5

#first round: test for differences among the clusters
#second round: test for differences among the clusters for the different conditions


print('Differential expression test on glutamatergic neurons.')
  
tmp_idx <- cellData$cell_type %in% "Neurons_Glutamatergic"
cellData_tmp <- cellData[tmp_idx,]
cellData_tmp$louvain <- factor(cellData_tmp$louvain)
sparse_mat_tmp <- sparse_mat[,tmp_idx]

#subcell types of interest are the louvain clusters that make up this cell type
cluster_OI <- levels(cellData_tmp$louvain)


#cell_idx <- cellData_tmp$louvain %in% cell_types
data_test <- as.matrix(sparse_mat_tmp)
#scale data to log-scran
data_test <- log(data_test %*% diag(1/cellData_tmp$size_factors) + 1)
rownames(data_test) <- genes

keep_genes <- rowMeans(data_test)>0.05
print(paste0('Keep ', sum(keep_genes), ' genes.'))
data_test <- data_test[keep_genes,]  

for (i in cluster_OI){
  cell_types <-i
  print(paste0('Test in cluster ', cell_types, '.'))
  
  condition_test <- cellData_tmp$louvain
  levels(condition_test) <- c(levels(condition_test), 'rest')
  condition_test[!condition_test == cell_types] <- 'rest'
  condition_test <- droplevels(condition_test)
  
  #test one cluster vs all_others
  rep_test <- cellData_tmp$sample
  rep_test <- droplevels(rep_test)
  #ct_test <- cellData$progen_test_type[cellData_tmp$progen_test_type %in% cell_types]
  #ct_test <- droplevels(ct_test)
  #cell_types_test <- levels(ct_test)
  #levels(ct_test) <- c('ct1', 'ct2')
  #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
  #?limma
  #prepare model matrix
  
  dge <- edgeR::DGEList(data_test,sample=rep_test, group=condition_test)
  design <- model.matrix(~condition_test + rep_test)
  colnames(design) <- gsub('condition_test','', colnames(design))
  #contrasts <- makeContrasts(ct1-ct2, levels = design)
  y <- new("EList")
  y$E <- dge
  rm(dge)
  
  fit <- limma::lmFit(y, design = design)
  rm(y)
  fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
  #fit2 <- contrasts.fit(fit, contrasts)
  rest <- limma::topTable(fit,coef= 2,
                          number = Inf,adjust.method='BH')
  
  subset_rest <- rest[abs(rest$logFC)>logFC_thresh & rest$adj.P.Val<p_adj_thresh,]
  title <- paste0(colnames(fit$coefficients)[2],'_vs_', 
                  names(table(condition_test))[!names(table(condition_test))==colnames(fit$coefficients)[2]], 
                  '_neuron_glu'
  )
  write.csv(x = subset_rest, 
            file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
  
  #create volcano plot
  data_res <- subset_rest %>% select(logFC, adj.P.Val)
  data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
  
  g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
    labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
    theme_classic()
  ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
         plot = g0, width = 6, height=6)
}

#######
# second round
#######
logFC_thresh <- 0.05
p_adj_thresh <- 1e-5

neuron_idx <- cellData$cell_type %in% "Neurons_Glutamatergic"
cellData_neuron <- cellData[neuron_idx,]
cellData_neuron$louvain <- factor(cellData_neuron$louvain)
sparse_mat_neuron <- sparse_mat[,neuron_idx]

#
ggplot(cellData_neuron, aes(louvain, n_counts, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cellData_neuron, aes(louvain, n_genes, fill=condition)) +geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#subcell types of interest are the louvain clusters that make up this cell type
cluster_OI <- levels(cellData_neuron$louvain)


#select by baseline

for (baseline_idx in levels(cellData_neuron$baseline)){
  print(paste0('Differential expression test on ', baseline_idx))
  tmp_idx <- cellData_neuron$baseline %in% baseline_idx
  cellData_tmp <- cellData_neuron[tmp_idx,]
  #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
  sparse_mat_tmp <- sparse_mat_neuron[,tmp_idx]
  
  for (i in cluster_OI){
    cell_types <-i
    print(paste0('Test in cluster ', cell_types, ' cells.'))
    
    cell_idx <- cellData_tmp$louvain %in% cell_types
    data_test <- as.matrix(sparse_mat_tmp[,cell_idx])
    #scale data to log-scran
    data_test <- log(data_test %*% diag(1/cellData_tmp$size_factors[cell_idx]) + 1)
    rownames(data_test) <- genes
    
    keep_genes <- rowMeans(data_test)>0.05
    print(paste0('Keep ', sum(keep_genes), ' genes.'))
    data_test <- data_test[keep_genes,]  
    
    condition_test <- cellData_tmp$adult_stress[cellData_tmp$louvain %in% cell_types]
    condition_test <- droplevels(condition_test)
    
    rep_test <- cellData_tmp$sample[cellData_tmp$louvain %in% cell_types]
    rep_test <- droplevels(rep_test)
    #ct_test <- cellData$progen_test_type[cellData_tmp$progen_test_type %in% cell_types]
    #ct_test <- droplevels(ct_test)
    #cell_types_test <- levels(ct_test)
    #levels(ct_test) <- c('ct1', 'ct2')
    #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
    #?limma
    #prepare model matrix
    
    dge <- edgeR::DGEList(data_test,sample=rep_test, group = condition_test)
   # design <- model.matrix(~condition_test + rep_test)
    design <- model.matrix(~condition_test)
    colnames(design) <- gsub('condition_test','', colnames(design))
    #contrasts <- makeContrasts(ct1-ct2, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef= 2,
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>logFC_thresh & rest$adj.P.Val<p_adj_thresh,]
    title <- paste0(baseline_idx, '_', colnames(fit$coefficients)[2],'_vs_', 
                    names(table(condition_test))[!names(table(condition_test))==colnames(fit$coefficients)[2]], 
                    '_', cell_types
    )
    write.csv(x = subset_rest, 
              file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
    
    #create volcano plot
    data_res <- subset_rest %>% select(logFC, adj.P.Val)
    data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
    
    g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
      labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
      theme_classic()
    ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
           plot = g0, width = 6, height=6)
  }
}

#select by adult_stress

for (stress_idx in levels(cellData$adult_stress)){
  print(paste0('Differential expression test on ', stress_idx))
  tmp_idx <- cellData_neuron$adult_stress %in% stress_idx
  cellData_tmp <- cellData_neuron[tmp_idx,]
  #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
  sparse_mat_tmp <- sparse_mat_neuron[,tmp_idx]
  
  for (i in cluster_OI){
    cell_types <-i
    print(paste0('Test in cluster ', cell_types, ' cells.'))
    
    cell_idx <- cellData_tmp$louvain %in% cell_types
    data_test <- as.matrix(sparse_mat_tmp[,cell_idx])
    #scale data to log-scran
    data_test <- log(data_test %*% diag(1/cellData_tmp$size_factors[cell_idx]) + 1)
    rownames(data_test) <- genes
    
    keep_genes <- rowMeans(data_test)>0.05
    print(paste0('Keep ', sum(keep_genes), ' genes.'))
    data_test <- data_test[keep_genes,]  
    
    condition_test <- cellData_tmp$baseline[cellData_tmp$louvain %in% cell_types]
    condition_test <- droplevels(condition_test)
    
    rep_test <- cellData_tmp$sample[cellData_tmp$louvain %in% cell_types]
    rep_test <- droplevels(rep_test)
    #ct_test <- cellData$progen_test_type[cellData_tmp$progen_test_type %in% cell_types]
    #ct_test <- droplevels(ct_test)
    #cell_types_test <- levels(ct_test)
    #levels(ct_test) <- c('ct1', 'ct2')
    #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
    #?limma
    #prepare model matrix
    
    dge <- edgeR::DGEList(data_test,sample=rep_test, group = condition_test)
    #design <- model.matrix(~condition_test + rep_test)
    design <- model.matrix(~condition_test)
    colnames(design) <- gsub('condition_test','', colnames(design))
    #contrasts <- makeContrasts(ct1-ct2, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef= 2,
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>logFC_thresh & rest$adj.P.Val<p_adj_thresh,]
    title <- paste0(stress_idx, '_', colnames(fit$coefficients)[2],'_vs_', 
                    names(table(condition_test))[!names(table(condition_test))==colnames(fit$coefficients)[2]], 
                    '_', cell_types
    )
    write.csv(x = subset_rest, 
              file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
    
    #create volcano plot
    data_res <- subset_rest %>% select(logFC, adj.P.Val)
    data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
    
    g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
      labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
      theme_classic()
    ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
           plot = g0, width = 6, height=6)
  }
}
