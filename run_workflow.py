import yaml
from analysis.preprocessing import load_gene_data, load_metadata
from analysis.processing import run_qc_norm
from analysis.plots import plot_qc_metrics

def main():
    """
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
    plot_qc_metrics(adata)

    #load in qc params
    n_feature_min = params['qc']['n_feature_min']
    n_feature_max = params['qc']['n_feature_max']
    n_count_min = params['qc']['n_count_min']
    n_count_max = params['qc']['n_count_max']
    percent_mt_max = params['qc']['percent_mt_max']

    #perform qc and normalization
    print("Performing QC and normalization...")
    adata, metadata = run_qc_norm(adata, metadata, n_feature_min, n_feature_max, n_count_min, n_count_max, percent_mt_max)

    print(adata)
    print(metadata)
    print(adata.obs)


if __name__ == '__main__':
    main() 