import scanpy as sc
import warnings

def run_qc_norm(adata, metadata, n_feature_min, n_feature_max, n_count_min, n_count_max, percent_mt_max):
    """
    """

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
    print('Number of TNBC samples pass QC:', len(adata[metadata['subtype'] == "TNBC"]))
    print('Number of HER2+ samples pass QC:', len(adata[metadata['subtype'] == "HER2+"]))
    print('Number of ER+ samples pass QC:', len(adata[metadata['subtype'] == "ER+"]))

    #perform normalization (CPM)
    warnings.simplefilter("ignore", UserWarning)
    sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True)
    metadata = adata.obs

    return adata, metadata