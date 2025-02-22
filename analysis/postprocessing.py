import scanpy as sc
from analysis.plots import plot_cell_labels

def umap_labels(adata, all_DEFs):
    """
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
        '0': 'Adipocytes',
        '1': 'T Cells',
        '2': 'B Cells',
        '3': 'Dendritic Cells',
        '4': 'Cardiomyocytes',
        '5': 'Fibroblasts',
        '6': 'Oligodendrocytes',
        '7': 'NK Cells',
        '8': 'Smooth Muscle Cells',
        '9': 'Basal Keratinocytes',
        '10': 'Epithelial Cells',
        '11': 'Suprabasal Keratinocytes',
        '12': 'Grandular Cells',
        '13': 'Enterocytes',
        '14': 'Kupffer Cells',
        '15': 'Plasma Cells',
        '16': 'Hofbauer Cells'
    }

    #replace metadata cell type labels with above labels
    adata.obs['cell_type'] = adata.obs['leiden'].map(cell_labels)

    #plot umap with all cell labels
    plot_cell_labels(adata)

    return adata

#cluster cell labels
#0: VWF - ADIPOCYTES OR ENDOTHELIAL CELLS
#1: GPR171 - T CELLS
#2: MS4A1 - B cells
#3: GZMB - DENDRITIC OR NK CELLS
#4: CRYAB -  CARDIOMYOCYTES OR PROXIMAL TUBULAR CELLS OR SKELETAL MYOCYTES
#5: CTSK - fibroblasts
#6: SCGB2B2/DKK1 - OLIGODENDROCYTES OR ENDOMETRIAL STROMAL CELLS
#7: IGKV1-5 - PLASMA CELLS
#8: NDUFA4L2 - SMOOTH MUSCLE CELLS OR KERATINOCYTES 
#9: KRT14 - BASAL KERATINOCYTES
#10: CD24 - EPITHELIAL CELLS
#11: S100A8 - SUPRABASAL KERATINOCYTES
#12: ANKRD30A - BREAST GLANDULAR CELLS
#13: AGR3 - PROXIMAL ENTEROCYTES
#14: HMOX1 - KUPFFER CELLS OR TROPHOBLASTS
#15: IGHG1 - PLASMA CELLS
#16: FTL - HOFBAUER CELLS