import pandas as pd
import scanpy as sc
import shutil
import os

def load_gene_data(raw_gene_data_path, raw_barcodes_path, raw_mtx_path, updated_data_path):
    """
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
    #df_genes_updated = pd.read_csv("./data/updated_genes.tsv", sep='\t', header=None)
    gene_data = sc.read_10x_mtx(updated_data_path, var_names='gene_symbols', cache=False)

    return gene_data

def load_metadata(raw_metadata_path, updated_data_path, gene_data):
    """
    """

    #read in metadata
    df_metadata = pd.read_csv(raw_metadata_path)

    #reformat metadata to match scanpy requirements
    new_headers = ['cell_id'] + list(df_metadata.columns[1:])
    df_metadata.columns = new_headers

    #save updated metadata
    df_metadata.to_csv(updated_data_path + 'metadata.csv', index=False)

    metadata = gene_data.obs.join(df_metadata.set_index('cell_id'))

    return metadata