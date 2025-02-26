import pandas as pd
import scanpy as sc
import shutil
import os
import warnings
from analysis.plots import plot_highvarfeats

def load_gene_data(raw_gene_data_path, raw_barcodes_path, raw_mtx_path, updated_data_path):
    """
    Loads in provided raw dataset to update sparse matrix, barcodes, and genes files into scanpy format

    Args:
        raw_gene_data_path (str): path of raw count_matrix_genes.tsv file
        raw_barcodes_path (str): path of raw count_matrix_barcodes.tsv file
        raw_mtx_path (str) : path of raw count_matrix_sparse.mtx file
        updated_data_path (str): path for updated genes, barcodes, and matrix files to be saved to

    Returns:
        adata (AnnData object): scanpy object in dataframe format containing gene expression values for each cell/sample
    """

    #read in gene data
    df_genes = pd.read_csv(raw_gene_data_path, sep='\t', header=None, names=['GeneName'])

    #reformat gene data to match scanpy requirements
    df_genes['GeneID'] = ['Gene_' + str(idx) for idx in range(len(df_genes))]
    df_genes = df_genes[['GeneID', 'GeneName']] 

    #save updated genes.tsv file
    df_genes.to_csv(updated_data_path + "genes.tsv", sep='\t', index=False, header=False)

    #copy and move matrix and barcode files
    shutil.copy(raw_barcodes_path, os.path.join(updated_data_path, 'barcodes.tsv'))
    shutil.copy(raw_mtx_path, os.path.join(updated_data_path, 'matrix.mtx'))

    #load and read in updated data files (genes, matrix, barcodes)
    adata = sc.read_10x_mtx(updated_data_path, var_names='gene_symbols', cache=True)

    return adata

def load_metadata(raw_metadata_path, updated_data_path, adata):
    """
    Loads in updated dataset to extract and append respective metadata

    Args:
        raw_metadata_path (str): path of raw metadata.csv file
        updated_data_path (str): path for updated genes, barcodes, and matrix files to be saved to
        adata (AnnData object): scanpy data object in dataframe format containing gene expression data per cell and metadata

    Returns:
        metadata (AnnData object): scanpy object in dataframe format containing only relevant metadata for each cell/sample
    """

    #read in metadata
    metadata = pd.read_csv(raw_metadata_path)

    #reformat metadata to match scanpy requirements
    new_headers = ['cell_id'] + list(metadata.columns[1:])
    metadata.columns = new_headers

    #save updated metadata
    metadata.to_csv(updated_data_path + 'metadata.csv', index=False)

    adata.obs = adata.obs.join(metadata.set_index('cell_id'))
    metadata = adata.obs

    #calculate % mitochondrial DNA per sample
    mito_genes = adata.var_names.str.startswith('MT-')
    metadata['percent_mt'] = adata[:, mito_genes].X.sum(axis=1)/adata.X.sum(axis=1)*100

    return metadata

def run_qc_norm(adata, metadata, n_feature_min, n_feature_max, n_count_min, n_count_max, percent_mt_max):
    """
    Loads in updated dataset and metadata to perform filtering and normalization

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing gene expression data per cell and metadata
        metadata (AnnData object): scanpy object in dataframe format containing only relevant metadata for each cell/sample
        n_feature_min (int): min threshold for feature/gene count of cell
        n_feature_max (int): max threshold for feature/gene count of cell
        n_count_min (int): min threshold for RNA count of cell
        n_count_max (int): max threshold for RNA count of cell
        percent_mt_max (int): max threshold for percent mitochondrial DNA to be present in the cell

    Returns:
        adata (AnnData object): scanpy data object in dataframe format containing filtered gene expression data per cell and metadata
        metadata (AnnData object): scanpy object in dataframe format containing only relevant filtered metadata for each cell/sample
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

    return adata, metadata

def select_feats(adata, n_HVFs, n_feats):
    """
    Loads in filtered and normalized data and perform highly variable feature extraction

    Args:
        adata (AnnData object): scanpy data object in dataframe format containing gene expression data per cell and metadata
        n_HVFs (int): number of features/genes to be included for analysis
        n_feats (int): number of most differnetially expressed features to be displayed/plotted

    Returns:
        adata (AnnData object): updated scanpy object in dataframe format containing only HVGs (default 2000)
    """

    #log transform of data to prep for PCA
    sc.pp.log1p(adata)

    #compute top 2000 highly variable features
    sc.pp.highly_variable_genes(adata, n_top_genes=n_HVFs, flavor='seurat', batch_key=None)

    #call fxn to display highly variable features
    plot_highvarfeats(adata)

    #subset data by highly variable features
    adata = adata[:, adata.var['highly_variable']]
    top_n_HVFs = adata.var.sort_values(by='dispersions_norm', ascending = False).head(n_feats)

    #display top genes/features
    print(f'Top {n_feats} highly variable features:')
    print(top_n_HVFs)

    #scale data
    print("Scaling data via linear transformation...")
    sc.pp.scale(adata)

    return adata