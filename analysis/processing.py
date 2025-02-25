import scanpy as sc
from analysis.plots import plot_pca, plot_umap
import pandas as pd
import warnings

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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

    #extract DEFs from rank_genes_groups
    DEF_result = adata.uns['rank_genes_groups']
    DEF_groups = DEF_result['names'].dtype.names
    all_DEFs = set()

    #get PC/group-specific markers and metrics
    for group in DEF_groups:
        group_DEFs = DEF_result['names'][group]
        group_scores = DEF_result['scores'][group]
        group_pvals = DEF_result['pvals'][group]
        all_DEFs.update(DEF_result['names'][group])

        #append and sort via pandas df
        DEF_df = pd.DataFrame({'Marker': group_DEFs, 'Score':group_scores, 'p-value': group_pvals})
        DEF_df_sorted = DEF_df.sort_values(by='Score', ascending=False)

        print(f'Top differentially expressed markers for PC {group}:') 
        print(DEF_df_sorted.head(n_genes))       

    #plot differentially expressed markers
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False)

    return adata, all_DEFs