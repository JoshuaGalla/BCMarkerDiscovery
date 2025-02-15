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
        ax = sc.pl.scatter(adata, x='nCount_RNA', y='nFeature_RNA', color = 'subtype', title = f'RNA Count by Gene Count for each Subtype', show = False)
        ax.set_xlabel('Total RNA Count')
        ax.set_ylabel('Total Gene Count')
        ax.text(0.50, 0.90, f'R^2 = {rsquared:.3f}', transform = ax.transAxes)
    plt.show()

    #RNA count by MT DNA percentage
    r, _ = pearsonr(adata.obs['nCount_RNA'], adata.obs['percent_mt'])
    rsquared = r**2
    with rc_context({"figure.figsize": (8, 6)}):
        ax = sc.pl.scatter(adata, x='nCount_RNA', y='percent_mt', color = 'subtype', title = f'RNA Count by percent_mt for each Subtype', show = False)
        ax.set_xlabel('Total RNA Count')
        ax.set_ylabel('Total percent MT DNA (%)')
        ax.text(0.50, 0.90, f'R^2 = {rsquared:.3f}', transform = ax.transAxes)
    plt.show()



