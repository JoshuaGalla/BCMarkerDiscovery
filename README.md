# **Read Me:**

Repository for BCMarkerDiscovery, an scRNA-seq analysis of breast cancer subtypes.

# **Project Overview:**

The goal of the project detailed in this repository is to perform a comprehensive analysis of the tumor microenvironment of three separate breast cancer subtypes - HER2+, ER+, and TNBC - from scRNA-seq data. This analysis consists of three stages defined by preprocessing (dataset loading, QC, normalization, highly variable feature selection, and scaling), processing (PCA, clustering, differentially expressed feature extraction, and UMAP plotting), and postprocessing (cell-specific UMAP labelling and plotting). The overall aim of this project is to identify cellular biomarkers that could provide potential diagnostic insight for differentiating between breast cancer subtypes. 

Additional information on the biological relevance of this analysis and in-depth step-by-step details of the workflow can be found in the *analysis_report.ipynb* notebook. All PC identifications were made by cross-referencing differentially expressed features with their single-cell identifications from [The Human Protein Atlas](https://www.proteinatlas.org/).

# **Data:**

The dataset used in this project can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) and consists of scRNA-seq data from 26 combined HER2+, ER+, and TNBC-specific patient samples. The supplementary file *GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz* can be found at the bottom of the GEO page and is the basis for this analysis.

The dataset should be downloaded and unzipped, and its contents placed into its own directory labeled *data/Wu_etal_2021_BRCA_scRNASeq* at the same file directory level as this *README.md* file and the *run_workflow.py* script. Within this new directory, ensure that the following files are present upon unzipping the downloaded dataset:
1) count_matrix_barcodes.tsv
2) count_matrix_genes.tsv
3) count_matrix_sparse.mtx
4) metadata.csv

Below is a visual example of what the directory structure should look like:
```
--BCMarkerDiscovery
    |--analysis
        |--...
    |--notebooks
        |--...
    |--parameters
        |--params.yaml
    |--data
        |--Wu_etal_2021_BRCA_scRNASeq
            |--count_matrix_barcodes.tsv
            |--count_matrix_genes.tsv
            |--count_matrix_sparse.mtx
            |--metadata.csv
    |--README.md
    |--requirements.txt
    |--run_workflow.py
```

# **Running the Workflow:**

To run the model, navigate to this directory and enter the following:

`python run_workflow.py`

Parameters for adjusting the raw/updated dataset file paths, QC thresholds, plot settings, number of PCs, and more can be adjusted in the _params.yaml_ file.

# **Output:**

Running the script *run_workflow.py* using the directions listed in the __Running the Workflow__ section above will display a series of outputs to the user's terminal and pop-up plots, including:

1) Number of samples that passed QC
2) Highly variable features across subtypes
3) PC elbow plot, PCA subplot for all subtypes, and differentially expressed markers amongst the top PCs
4) Most differentially expressed markers across all selected PCs
5) Marker-specific UMAP
6) Overall UMAP with cell-specific identifications

# **Dependencies:**

Prerequisites for running this model are included in the _requirements.txt_ file.

# **License:**

This project is licensed under the MIT license. See the LICENSE file/tab for more info.

# **Acknowledgements:**

This project was inspired by the work of [Ido Nofech-Mozes et al.](https://www.nature.com/articles/s41467-023-37353-8) and the [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) developed by Rahul Satija, Satija Lab, and Collaborators.





