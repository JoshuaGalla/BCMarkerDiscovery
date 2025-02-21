import scanpy as sc
from analysis.plots import plot_pca, plot_umap
import numpy as np

def run_pca(adata, n_vars, pc_range):
    """
    """

    #applying PCA
    sc.tl.pca(adata, svd_solver='auto')

    #plot top PCs, dimensionality heatmap, elbow plot
    plot_pca(adata, n_vars, pc_range)

    return adata

def cluster_data(adata, n_neighbors, n_pcs, resolution, min_dist, spread, display_UMAP_unlabeled):
    """
    """

    #compute nearest neighbor distances and matrix using Gaussian kernel for KNN
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    #cluster cells into subgroups via leiden algorithm
    sc.tl.leiden(adata, n_iterations=2, resolution=resolution, flavor='igraph', directed=False)

    #embed neighborhood graph using UMAP
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)

    #plot UMAP
    plot_umap(adata, n_pcs, display_UMAP_unlabeled)

    return adata

def find_DEFs(adata, n_genes):
    """
    """

    #find differentially expressed markers between PCs
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

    #plot differentially expressed markers
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False)

    #plot specific marker of interest on UMAP
    sc.pl.umap(adata, color='VWF', use_raw=False)
    sc.pl.umap(adata, color='GPR171', use_raw=False)
    sc.pl.umap(adata, color='GZMB', use_raw=False)
    sc.pl.umap(adata, color='MS4A1', use_raw=False)
    sc.pl.umap(adata, color='CRYAB', use_raw=False)
    sc.pl.umap(adata, color='SCGB1B2P', use_raw=False)
    sc.pl.umap(adata, color='IGKV1-5', use_raw=False)
    sc.pl.umap(adata, color='NDUFA4L2', use_raw=False)
    sc.pl.umap(adata, color='KRT14', use_raw=False)
    sc.pl.umap(adata, color='CD24', use_raw=False)
    sc.pl.umap(adata, color='S100A8', use_raw=False)
    sc.pl.umap(adata, color='ANKRD30A', use_raw=False)
    sc.pl.umap(adata, color='AGR3', use_raw=False)
    sc.pl.umap(adata, color='HMOX1', use_raw=False)
    sc.pl.umap(adata, color='IGHG1', use_raw=False)
    sc.pl.umap(adata, color='FTL', use_raw=False)

    return adata