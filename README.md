# **Read Me:**

Repository for BIOINF 545 final project.

# **Project Overview:**

The goal of this project is to perform a comprehensive analysis of the tumor microenvironment of three separate breast cancer subtypes - HER2+, ER+, and TNBC - from scRNA-seq data. This analysis consists of data QC and normalization, variable feature selection, scaling and linear dimensionality reduction (PCA), clustering and non-linear dimensionality reduction (UMAP), and cell-type population labeling using R's Seurat Object. 

This project was inspired by the work of [Ido Nofech-Mozes et al.](https://www.nature.com/articles/s41467-023-37353-8), in which they developed a new model for classifying cell types in the tumor microenvironment across cancer types called scATOMIC, which improves cellular classification and analysis of the TME setting. Previously, cell-type classifications within the TME have been difficult due to high heterogeneity among the same tissue type and low transcriptomic variation among specialized immune cells. scATOMIC aims to bridge this knowledge gap and its pipeline has been specifically extended/built upon for breast cancer classification.

# **Data:**

This project specifically builds upon the work above by attempting to replicate their results using a [separate dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) of 26 combined HER2+, ER+, and TNBC-specific patient samples. The supplementary file `GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz` can be found at the bottom of the GEO page and is the basis for this analysis.

The dataset should be downloaded and unzipped, and its contents placed into its own directory labeled `data` at the same level as `analysis` and `report`. Within this new directory, the following files should be present:
1) barcodes.tsv
2) genes.tsc
3) matrix.mtx
4) metadata.csv

# **Workflow:**

An R script for conducting this project's analysis is described in the `analysis/project_analysis.R` file. A subsequent markdown file and formal report that contains this raw code in combination with further annotations, links to literature supporting findings, figures/plots, and continued discussion and future directions can be seen in the files `report/project_analysis_markdown.Rmd` and `report/project_analysis_report.html` respectively. 

Due to the file size viewing limitations on GitHub, the formal report file `report/project_analysis_report.html` may need to be downloaded in order to be viewed. 

A Python alternative to the analysis performed in `analysis/project_analysis.R` is currently being developed, and is viewable in the Jupyter notebook `analysis/project_analysis.ipynb`. In the near future, a .py script will also be implemented to run this analysis and display plots of interest for a specified BC subtype.





