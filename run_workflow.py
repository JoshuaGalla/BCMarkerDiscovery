import yaml
from analysis.preprocessing import load_gene_data, load_metadata

def main():
    """
    """

    #read in parameters
    print("Reading dataset...")
    with open('parameters/params.yaml', 'r') as file:
        params = yaml.safe_load(file)

    #set raw and updated data paths
    raw_gene_data_path = params['raw_data']['raw_gene_data_path']
    raw_barcodes_path = params['raw_data']['raw_barcodes_path']
    raw_mtx_path = params['raw_data']['raw_mtx_path']
    updated_data_path = params['updated_data']['updated_data_path']

    #load in gene data, metadata
    gene_data = load_gene_data(raw_gene_data_path, raw_barcodes_path, raw_mtx_path, updated_data_path)

    #load in metadata
    raw_metadata_path = params['raw_data']['raw_metadata_path']
    metadata = load_metadata(raw_metadata_path, updated_data_path, gene_data)

    return gene_data, metadata

if __name__ == '__main__':
    main() 