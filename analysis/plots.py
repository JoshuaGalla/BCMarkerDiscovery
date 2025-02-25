import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scipy.stats import pearsonr

def plot_qc_metrics(adata):
    """
    """

    #gene count distribution
    with rc_context({"figure.figsize": (8,6)}):
        ax = sc.pl.violin(adata, ['nFeature_RNA'], groupby='subtype', inner='box', stripplot=True, show=False)
        ax.set_title('Gene Count Distribution by Subtype')
        ax.set_xlabel('Subtype')
        ax.set_ylabel('Gene Count')
    plt.show()

    #RNA count distribution
    with rc_context({"figure.figsize": (8,6)}):
        ax = sc.pl.violin(adata, ['nCount_RNA'], groupby='subtype', inner='box', stripplot=True, show=False)
        ax.set_title('RNA Count Distribution by Subtype')
        ax.set_xlabel('Subtype')
        ax.set_ylabel('RNA Count')
    plt.show()

    #percent MT DNA distribution
    with rc_context({"figure.figsize": (8,6)}):
        ax = sc.pl.violin(adata, ['percent_mt'], groupby='subtype', inner='box', stripplot=True, show=False)
        ax.set_title('Percent Mitochondrial DNA Distribution by Subtype')
        ax.set_xlabel('Subtype')
        ax.set_ylabel('Percent MT DNA (%)')
    plt.show()

    #RNA count by gene count
    r, _ = pearsonr(adata.obs['nCount_RNA'], adata.obs['nFeature_RNA'])
    rsquared = r**2
    with rc_context({"figure.figsize": (8, 6)}):
        ax = sc.pl.scatter(adata, x='nCount_RNA', y='nFeature_RNA', color = 'subtype', title = f'RNA Count by Gene Count per Subtype', show = False)
        ax.set_xlabel('Total RNA Count')
        ax.set_ylabel('Total Gene Count')
        ax.text(0.50, 0.90, f'R^2 = {rsquared:.3f}', transform = ax.transAxes)
    plt.show()

    #RNA count by MT DNA percentage
    r, _ = pearsonr(adata.obs['nCount_RNA'], adata.obs['percent_mt'])
    rsquared = r**2
    with rc_context({"figure.figsize": (8, 6)}):
        ax = sc.pl.scatter(adata, x='nCount_RNA', y='percent_mt', color = 'subtype', title = f'RNA Count by MT DNA per Subtype', show = False)
        ax.set_xlabel('Total RNA Count')
        ax.set_ylabel('Total percent MT DNA (%)')
        ax.text(0.50, 0.90, f'R^2 = {rsquared:.3f}', transform = ax.transAxes)
    plt.show()

def plot_highvarfeats(adata):
    """
    """

    #display top 2000 highly variable features
    with rc_context({"figure.figsize": (8,6)}):
        sc.pl.highly_variable_genes(adata, show=True)

def plot_pca(adata, n_vars, pc_range):
    """
    """

    #display elbow plot to determine dimensionality
    with rc_context({"figure.figsize": (8,6)}):
        sc.pl.pca_variance_ratio(adata, show=True)

    #display PCA subplot
    with rc_context({"figure.figsize": (8,6)}):
        ax = sc.pl.pca(adata, color='subtype', show=False)
        ax.set_title('PCA for Highest PCs Across All Subtypes')
    plt.show()

    #display top n PCs
    sc.pl.pca_loadings(adata, components=pc_range, n_points=n_vars)

def plot_umap(adata, n_pcs, display_UMAP_unlabeled):
    """
    """

    #display UMAP of embedded clusters for n PCs (labeled by subtype)
    with rc_context({"figure.figsize": (8,6)}):
        ax=sc.pl.umap(adata, color=['subtype'], show=False)
        ax.set_title(f'UMAP of {n_pcs} PC Clusters Labeled by Subtype')
    plt.show()

    if display_UMAP_unlabeled == True:
        #display UMAP of embedded clusters for n PCs (not labeled)
        with rc_context({"figure.figsize": (8,6)}):
            ax=sc.pl.umap(adata, color='leiden', show=False)
            ax.set_title(f'UMAP of {n_pcs} PC Clusters Unlabeled')
        plt.show()

def plot_cell_labels(adata):
    """
    """

    with rc_context({"figure.figsize": (8,6)}):
        sc.pl.umap(adata, color='cell_type', title = 'UMAP of Labeled PC Clusters', legend_loc = 'on data', legend_fontsize=6, legend_fontoutline=2)

   
    




