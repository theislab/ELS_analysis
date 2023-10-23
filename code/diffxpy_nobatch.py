import os 
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

import scanpy as sc
import anndata as ann
import diffxpy.api as de
import batchglm.api as glm
import scanpy.external as sce

from scipy.sparse import csr_matrix

sc.settings.verbosity = 3
sc.logging.print_versions()
print(de.__version__)

adata = sc.read('../data/processed/adata_annotated.h5ad')
counts = sc.read('../data/processed/data_norm.h5ad')
counts.obs = adata.obs.copy()
gene_idx = (counts.var.index.isin(adata.var.index))

de_results_celltype = []
de_results_louvain = []
for test in ['BL', 'ELS']:
    adata_loop = adata[adata.obs['baseline'] == test]
    counts_loop = counts[counts.obs['baseline'] == test]
    de_results = {'test': test, 'pair': 'baseline'}
    de_results_coarse = {'test': test, 'pair': 'baseline'}
                
    for clust in adata_loop.obs['louvain'].cat.categories:
        adata_tmp = adata_loop[adata_loop.obs['louvain'] == clust,:].copy()
        counts_tmp = counts_loop[counts_loop.obs['louvain'] == clust,:].copy()
        print(f'In cluster {clust}:')
        # Filter out genes to reduce multiple testing burden
        sc.pp.filter_genes(counts_tmp, min_cells=np.ceil(0.05*counts_tmp.shape[0]))
        print(f'Testing {counts_tmp.n_vars} genes...')
        print("")
        try:
            coefs = de.utils.preview_coef_names(
                sample_description=counts_tmp.obs,
                formula="~1+condition", 
            )
            print(coefs)
            test_tmp = de.test.wald(
                data=counts_tmp.layers['counts'],
                formula_loc="~ 1 + condition",
                size_factors='size_factors',
                coef_to_test=[coefs[-1]],
                sample_description=counts_tmp.obs,
                gene_names=counts_tmp.var_names,
                noise_model='nb',
                dtype="float64",
            )
            de_results[clust] = test_tmp
            de_results['coef_test'] = coefs[-1]
        except (np.linalg.LinAlgError, ValueError) as e:
            print('+++++/// Not enough sample ids in cluster! ///+++++')
            pass

        #Store the results
        
    de_results_louvain.append(de_results)

##########################################################
for test in ['Stress', 'Ctrl']:
    adata_loop = adata[adata.obs['adult_stress'] == test]
    counts_loop = counts[counts.obs['adult_stress'] == test]
    de_results = {'test': test, 'pair': 'adult_stress'}

    for clust in adata_loop.obs['louvain'].cat.categories:
        adata_tmp = adata_loop[adata_loop.obs['louvain'] == clust,:].copy()
        counts_tmp = counts_loop[counts_loop.obs['louvain'] == clust,:].copy()
        print(f'In cluster {clust}:')
        # Filter out genes to reduce multiple testing burden
        
        sc.pp.filter_genes(counts_tmp, min_cells=np.ceil(0.05*counts_tmp.shape[0]))
        print(f'Testing {counts_tmp.n_vars} genes...')
        print("")
        try:
            coefs = de.utils.preview_coef_names(
                sample_description=counts_tmp.obs,
                formula="~1+baseline", 
            )
            print(coefs)
            test_tmp = de.test.wald(
                data=counts_tmp.layers['counts'],
                formula_loc="~ 1 + baseline",
                size_factors='size_factors',
                coef_to_test=[coefs[-1]],
                sample_description=counts_tmp.obs,
                gene_names=counts_tmp.var_names,
                noise_model='nb',
                dtype="float64"
            )
            de_results[clust] = test_tmp
            de_results['coef_test'] = coefs[-1]
        except (np.linalg.LinAlgError, ValueError) as e:
            print('+++++/// Not enough sample ids in cluster! ///+++++')
            pass
        #Store the results
        
    de_results_louvain.append(de_results)
  
    with open('../results/diffxpy_louvain_correct_noconstraints_coeftotest_2021.pickle', 'wb') as handle:
        pickle.dump(de_results_louvain, handle, protocol=pickle.HIGHEST_PROTOCOL)