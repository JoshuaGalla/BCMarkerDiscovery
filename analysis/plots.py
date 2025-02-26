import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scipy.stats import pearsonr

def plot_qc_metrics(adata):
    """
    Displays QC info required for informed preprocessing, including gene, RNA, and percent mitochondrial DNA counts per subtype

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing subtype label, gene/RNA/MT% counts, and other metadata

    Returns:
        None - plots QC metric graphs per subtype
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
    Displays normalized and non-normalized top 2000 most highly variable features across all three subtypes

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing gene expression values for each cell/sample across subtypes

    Returns:
        None - plots most highly variable features across all subtypes
    """

    #display top 2000 highly variable features
    with rc_context({"figure.figsize": (8,6)}):
        sc.pl.highly_variable_genes(adata, show=True)

def plot_pca(adata, n_vars, pc_range):
    """
    Displays data dimensionality info required for UMAP, including PC elbow plot, PCA clustered by subtype, and top PCs contributing to variance

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing gene expression values per cell, PC identification, and other relevant metadata
        n_vars (int): number of genes/feature to be included per PC
        pc_range (list): number/range of PCs with variable features to be displayed

    Returns:
        None - plots elbow plot, PCA subplot, and differential features per displayed PC
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
    Displays single-cell, subtype-specific, and overall UMAP of top 2000 HVGs

    Args:
        adata (AnnData object): scanpy data object in dataframe format
        n_pcs (int): number of PCs considered for clustering; determined visually by looking at elbow plot
        display_UMAP_unlabeled (boolean): indicates whether unlababeled raw UMAP is displayed for user to view

    Returns:
        None - plots subtype-specific and unlabeled UMAPs
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
    Displays UMAP with previously identiifed clusters labeled with cell type based on canonical marker expression
    
    Args:
        adata (AnnData object): scanpy data object in dataframe format
        
    Returns:
        None - plots UMAP containing clusters labeled by significant cell type markers
    """

    with rc_context({"figure.figsize": (8,6)}):
        sc.pl.umap(adata, color='cell_type', title = 'UMAP of Labeled PC Clusters', legend_loc = 'on data', legend_fontsize=6, legend_fontoutline=2)

   
    




