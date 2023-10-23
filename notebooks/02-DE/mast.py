#!/usr/bin/env python3

import os
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

import rpy2.rinterface_lib.callbacks
import logging
import warnings
from rpy2.robjects import pandas2ri
import anndata2ri

rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()

from rpy2 import robjects
def mast_de_groups(adata, groupby, save):
    '''Compute differential expression with the MAST package by treatment covariate within clusters provided as "groupby" and export as excel file'''
    
    #if groupby not in ['louvain_final', 'louvain_three']:
    #    raise NotImplementedError("This function is only implemented to be used with 'louvain_final' and 'louvain_three' group labels")
    
    # Define R function to run MAST analysis
    robjects.r('''
        mast_de_r <- function(data_mat, clusters, obs, var, groupby){
            print('Deploying to R...')
            #Prepare data sets for SingleCellExperiment data structure conversion
            obs['wellKey'] = row.names(obs)
            var['primerid'] = row.names(var)

            #Convert to SingleCellExperiment type
            sca <- FromMatrix(exprsArray=data_mat, cData=obs, fData=var)

            #Compute Gene detection rate
            cdr <- colSums(assay(sca)>0)
            colData(sca)$ngeneson = scale(cdr)

            colData(sca)$n_genes = scale(colData(sca)$n_genes)
            #Create a vector that will hold all the DE results
            output <- vector("list", length(clusters))

            count <- 0
            print('Begin computation...')
            #Loop over all louvain clusters
            for (i in clusters){
                count <- count+1
                print(i)
                #Create data subsets which should be used for testing
                if (groupby=='louvain') {
                    sca_sub <- sca[,colData(sca)$louvain==i]
                } else if (groupby=='louvain_r1') {
                    sca_sub <- subset(sca, with(colData(sca), louvain_r1==i))
                } else {
                    stop()
                }
                #Filter out non-expressed genes in the subset
                sca_sub <- sca_sub[rowSums(assay(sca_sub)) != 0, ]
                #Define & run hurdle model
                zlmCond <- zlm(formula = ~condition + n_genes + sample, sca=sca_sub)
                summaryCond <- summary(zlmCond, doLRT='conditionStress')
                summaryDt <- summaryCond$datatable
                result <- merge(summaryDt[contrast=='conditionStress' & component=='H',.(primerid, `Pr(>Chisq)`)], #p-vals
                                 summaryDt[contrast=='conditionStress' & component=='logFC', .(primerid, coef)], #logFC coefficients
                                 by='primerid') 
                #Correct for multiple testing (FDR correction) and filtering
                result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
                result[,coef:=result[,coef]/log(2)]
                names(result) <- c("gene", "pval", "log2FC", "qval")
                result = result[order(result$qval),]
                print('debug4')
                output[[count]] <- result
                print('debug5')
            }
            return(output)
        }
    ''')
    
    mast_de = robjects.globalenv['mast_de_r']
    

    
    #Create new Anndata object for use in MAST with non-batch corrected data as before
    adata_test = adata.copy()
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) 
    obs = adata_test.obs
    var = adata_test.var
    data_mat = adata_test.X.T
    clusters = list(adata_test.obs[groupby].cat.categories)
    
    expr_dict = {adata_test.var.index[i]:{} for i in range(adata_test.shape[1])}
    expr_dict_stress = {adata_test.var.index[i]:{} for i in range(adata_test.shape[1])}
    expr_dict_ctrl = {adata_test.var.index[i]:{} for i in range(adata_test.shape[1])}
    for clust in adata_test.obs[groupby].cat.categories:
        expr = np.mean(adata_test[adata_test.obs[groupby] == clust].X, axis=0)
        expr_stress = np.mean(adata_test[(adata_test.obs['condition']=='Stress') & (adata_test.obs[groupby] == clust)].X, axis=0)
        expr_ctrl = np.mean(adata_test[(adata_test.obs['condition']=='Control') & (adata_test.obs[groupby] == clust)].X, axis=0)
        for i, gene in enumerate(adata_test.var.index):
            expr_dict[gene][clust] = expr[i]
            expr_dict_stress[gene][clust] = expr_stress[i]
            expr_dict_ctrl[gene][clust] = expr_ctrl[i]
    
    result = mast_de(data_mat, clusters, obs, var, groupby)
    result = {clusters[i]:datframe for i, datframe in enumerate(result)}
    
    writer = pd.ExcelWriter(save, engine='xlsxwriter')
    print('Number of significant DE genes:')    
    for clust in clusters:
        result[clust]['meanExpr'] = [expr_dict[gene][clust] for gene in result[clust]['gene'].values]
        result[clust]['meanExprStress'] = [expr_dict_stress[gene][clust] for gene in result[clust]['gene'].values]
        result[clust]['meanExprCtrl'] = [expr_dict_ctrl[gene][clust] for gene in result[clust]['gene'].values]
        result[clust].to_excel(writer,sheet_name=str(clust))
        print(clust+':', np.sum([result[clust]['qval']<0.05]))

    writer.save()

    return result

adata = sc.read('../data/processed/adata_annotated.h5ad')
adata_raw = sc.read('../data/processed/data_norm.h5ad')
adata_raw.obs = adata.obs
sc.pp.highly_variable_genes(adata_raw, n_top_genes=4000, flavor='cell_ranger', subset=True)
adata_bl = adata_raw[adata_raw.obs['baseline']=='BL'].copy()
adata_bl.obs = adata_bl.obs.rename(columns={'sample_id':'sample'})
adata_bl.obs['condition'] = adata_bl.obs['condition'].cat.rename_categories({'BL_Ctrl':'Ctrl',
                                                                            'BL_Stress':'Stress'})

de = mast_de_groups(adata_bl, groupby='louvain', save='../results/mast_de_baseline.xlsx')