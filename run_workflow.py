import yaml
from analysis.preprocessing import load_gene_data, load_metadata, run_qc_norm, select_feats
from analysis.processing import run_pca, cluster_data, find_DEFs
from analysis.postprocessing import umap_labels
from analysis.plots import plot_qc_metrics

def main():
    """
    Runs scanpy scRNA-seq workflow to perform preprocessing, processing, and postprocessing of inputted dataset 

    Args:
        None

    Returns:
        None - displays analysis progress, metrics and highly variable markers, and relevant plots
    """

    #read in parameters
    with open('parameters/params.yaml', 'r') as file:
        params = yaml.safe_load(file)

    #set raw and updated data paths
    raw_gene_data_path = params['raw_data']['raw_gene_data_path']
    raw_barcodes_path = params['raw_data']['raw_barcodes_path']
    raw_mtx_path = params['raw_data']['raw_mtx_path']
    updated_data_path = params['updated_data']['updated_data_path']

    #load in gene data
    print("Reading sample data...")
    adata = load_gene_data(raw_gene_data_path, raw_barcodes_path, raw_mtx_path, updated_data_path)

    #load in metadata
    raw_metadata_path = params['raw_data']['raw_metadata_path']
    print("Reading metadata...")
    metadata = load_metadata(raw_metadata_path, updated_data_path, adata)

    #plot subtype qc metrics
    display_qc_metrics = params['plots']['display_qc_metrics']
    if display_qc_metrics == True:
        plot_qc_metrics(adata)

    #load in qc params
    n_feature_min = params['QC']['n_feature_min']
    n_feature_max = params['QC']['n_feature_max']
    n_count_min = params['QC']['n_count_min']
    n_count_max = params['QC']['n_count_max']
    percent_mt_max = params['QC']['percent_mt_max']

    #perform qc and normalization
    print("Filtering samples per QC metrics...")
    adata, metadata = run_qc_norm(adata, metadata, n_feature_min, n_feature_max, n_count_min, n_count_max, percent_mt_max)

    #identify (and display) highly variable features
    n_HVFs = params['features']['n_HVFs']
    n_feats = params['features']['n_feats']
    print("Identifying highly variable features...")
    adata = select_feats(adata, n_HVFs, n_feats)

    #run PCA
    n_vars = params['PCA']['n_vars']
    pc_range = params['PCA']['pc_range']
    print("Performing linear dimensionality reduction (PCA)...")
    adata = run_pca(adata, n_vars, pc_range)

    #cluster highly variable features
    n_neighbors = params['clustering']['n_neighbors']
    n_pcs = params['clustering']['n_pcs']
    resolution = params['clustering']['resolution']
    min_dist = params['UMAP']['min_dist']
    spread = params['UMAP']['spread']
    display_UMAP_unlabeled = params['plots']['display_UMAP_unlabeled']
    print('Clustering data and creating UMAP...')
    adata = cluster_data(adata, n_neighbors, n_pcs, resolution, min_dist, spread, display_UMAP_unlabeled)

    #find differentially expressed features
    n_genes = params['DEF']['n_genes']
    print("Finding differentially expressed features...")
    adata, all_DEFs = find_DEFs(adata, n_genes)

    #plot marker-specific umaps and label umap with cell types
    adata = umap_labels(adata, all_DEFs)

if __name__ == '__main__':
    main() 