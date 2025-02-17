import scanpy as sc
import pandas as pd
import warnings
from analysis.plots import plot_highvarfeats, plot_pca

def run_qc_norm(adata, metadata, n_feature_min, n_feature_max, n_count_min, n_count_max, percent_mt_max):
    """
    """

    #calculating total dataset size
    num_her = len(adata[metadata['subtype'] == "HER2+"])
    num_er = len(adata[metadata['subtype'] == "ER+"])
    num_tnbc = len(adata[metadata['subtype'] == "TNBC"])

    #filter samples based on QC metrics set in parameters
    adata = adata[(
        (metadata['nFeature_RNA']>n_feature_min) &
        (metadata['nFeature_RNA']<n_feature_max) &
        (metadata['nCount_RNA']>n_count_min) &
        (metadata['nCount_RNA']<n_count_max) &
        (metadata['percent_mt']<percent_mt_max)
    )]
    metadata = adata.obs

    #display qc metrics
    num_her_qc = len(adata[metadata['subtype'] == "HER2+"])
    num_er_qc = len(adata[metadata['subtype'] == "ER+"])
    num_tnbc_qc = len(adata[metadata['subtype'] == "TNBC"])
    print('Number of TNBC samples pass QC (total, %):', num_tnbc_qc, round((num_tnbc_qc/num_tnbc)*100, 2))
    print('Number of HER2+ samples pass QC (total, %):', num_her_qc, round((num_her_qc/num_her)*100, 2))
    print('Number of ER+ samples pass QC (total, %):', num_er_qc, round((num_er_qc/num_er)*100, 2))

    #perform normalization (CPM)
    warnings.simplefilter("ignore", UserWarning)
    print("Normalizing data by CPM...")
    sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True)
    #metadata = adata.obs

    return adata, metadata

def select_feats(adata, n_highvarfeats, n_feats):
    """
    """

    #log transform of data to prep for PCA
    sc.pp.log1p(adata)

    #compute top 2000 highly variable features
    sc.pp.highly_variable_genes(adata, n_top_genes=n_highvarfeats, flavor='seurat', batch_key=None)

    #call fxn to display highly variable features
    plot_highvarfeats(adata)

    #subset data by highly variable features
    adata = adata[:, adata.var['highly_variable']]
    top_n_feats = adata.var.sort_values(by='dispersions_norm', ascending = False).head(n_feats)

    #display top genes/features
    print(f'Top {n_feats} highly variable features:')
    print(top_n_feats)

    #scale data
    print("Scaling data via linear transformation...")
    sc.pp.scale(adata)
    print(adata.shape)

    return adata

def run_pca(adata, n_vars, pc_range):
    """
    """

    #applying PCA
    sc.tl.pca(adata, svd_solver='auto')

    #plot top PCs, dimensionality heatmap, elbow plot
    plot_pca(adata, n_vars, pc_range)

    return adata