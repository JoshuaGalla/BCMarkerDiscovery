import scanpy as sc
from analysis.plots import plot_cell_labels

def umap_labels(adata, all_DEFs):
    """
    labels UMAP clusters with cell type, depending on presence of canonical markers. Enable marker-specific search and UMAP display

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing gene expression data per cell and metadata
        all_DEFs (set): set of all differnetially expressed features that exist across top PCs; enables search function

    Returns:
        adata (AnnData object): updated scanpy object in dataframe format
    """

    #search and plot specific marker of interest on UMAP
    while True:
        DEF_name = input('Enter biomarker name to plot (or type "exit" to continue with analysis):')

        if DEF_name.lower() == 'exit':
            break

        #check if gene/input exists
        if DEF_name in all_DEFs:
            sc.pl.umap(adata, color=DEF_name, use_raw=False)

        else:
            print(f"Biomarker {DEF_name} not present in differentially expressed markers. Please try again")

    #assign cluster labels of specific cell types from DEF analysis
    cell_labels = {
        '0': 'Adipocytes', #VWF
        '1': 'T Cells', #GPR171
        '2': 'B Cells', #MS4A1
        '3': 'Dendritic Cells', #GZMB
        '4': 'Cardiomyocytes', #CRYAB
        '5': 'Fibroblasts', #CTSK
        '6': 'Oligodendrocytes', #SCGB2B2
        '7': 'NK Cells', #IGKV1-5
        '8': 'Smooth Muscle Cells', #NDUFA4L2
        '9': 'Basal Keratinocytes', #KRT14
        '10': 'Epithelial Cells', #CD24
        '11': 'Suprabasal Keratinocytes', #S100A8
        '12': 'Grandular Cells', #ANKRD30A
        '13': 'Enterocytes', #AGR3
        '14': 'Kupffer Cells', #HMOX1
        '15': 'Plasma Cells', #IGHG1
        '16': 'Hofbauer Cells' #FTL
    }

    #replace metadata cell type labels with above labels
    adata.obs['cell_type'] = adata.obs['leiden'].map(cell_labels)

    #plot umap with all cell labels
    plot_cell_labels(adata)

    return adata