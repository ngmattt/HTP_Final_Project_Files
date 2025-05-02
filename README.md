Title: HTP_Project

Author: Matthew Ng

Date: 4/30/25

Language: R 4.4.2, Python 3.11
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Dependencies for R: DESeq2, stringr, ggplot2, dplyr, clusterProfiler, org.Hs.eg.db, kibble, AnnotationDbi, ggrepel

Dependencies for Python: pandas, dumpy, matplotlib, scikit-learn (sklearn)

Package instructions for R: Packages were installed in RStudio using the following command:

install.packages("name of the package") <- to install the package

library(name of the package) <- to use the package

Package instructions for Python: Packages were installed in Pycharm using the following command:
pip install "package name" <- to install the package

import "package name" <- to use the package

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Input: "GSE169687_raw_counts_GRCh38.p13_NCBI.tsv.gz", "GSE169687_series_matrix(1).txt"; downloaded from GEO Database, "HTP_Project_ScriptR.R", "HTP_Project_Script_3.ipynb"

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Output: "ML_input_all_volcano_genes.csv" is a file output from "HTP_Project_ScriptR.R" that was used for the "HTP_Project_Script_3.ipynb". 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Description: "HTP_Project_ScriptR.R" is the R script that is used to process the raw data files that were downloaded from GEO database. The purpose of this script is to prepare the raw counts file and series matrix file for differential gene expression analysis and to output "ML_input_all_volcano_genes.csv" which contains the list of significant DEGs from the analysis. The output file is then used for the "HTP_Project_Script_3.ipynb" which is a Jupyter notebook file that contains the code for the machine learning portion of the project. Running the script will provide performance metrics for each model trained and tested. 
