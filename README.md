## "Immune checkpoint inhibitor response in sarcomas associates with immune infiltrates and increased expression of transposable elements and viral response pathways" 

# Data description

The data directory includes the source data and code necessary to recreate all main and supplemental figures. The data directory has two subdirectories: 'study' and 'TCGA'.
The 'study' subdirectory features analysis of the clinical dataset collected at MSKCC, while the 'TCGA' subdirectory represents the data set and analysis that replicate our observations.

All necessary data files to reproduce all the figures are available in the 'data/study' and 'data/TCGA' folders and include:



C_GAS_gene_list_KEGG

Clinical_phenotypes

Epigentic_pathway_genes

Gene_lists_for_ddGSEA.txt

GLMnet_Functions.R

NanostringCentroids_from_Petiprez_publication

Pathway_annotations_for_24_pathways.txt

TE_family_annotation


NOTE: study and TCGA folders do not include gene and Transposable element count matrices ("gene_expression counts" and "INTERGENIC_TE_normalized_expression_1052_repeats"
becaus of their size), and these files could be downloaded from Zenodo: https://zenodo.org/uploads/10313854 for TCGA data, and study data are avable with controlled access from 
dbgap acccession phs003284. 

# Code description 
All the costume code used to process RNA-seq data after intial processing by REDISCOVERTE can be found in the 'data/study/Scripts/' folder or 'data/TCGA/Scripts/' folder.

## Script_1.R 

This script performs a comprehensive analysis of our clinical raw RNA sequencing counts obtained from the REDISCOVERTE pipeline run. The script first normalizes and filters the counts, and then performs immune deconvolution using the MCP counter. The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major clusters: immune hot and cold.

The script generates several figures, including a heatmap of the immune clusters and relevant clinical correlates (Figure 1A,B,C), a Kaplan-Meier analysis and plot for Immune hot/cold groups and PFS (months) (Figure 1D), factor miner clustering of deconvoluted data (Supplemental Figure 1), a Cox Hazardous model analysis and plot (Supplemental Figure 2), and a matching of our immune hot and cold clusters with previously published Petitprez et.al clustering of immune cells (Supplemental Figure 3). Additionally, the script generates a heatmap of normalized intergenic TE counts with clinical correlates (Supplemental Figure 5).

The script also includes several supplementary figures, including a correlation between IKZF1 and B-cells (Supplemental Figure 6A).


## Script_2.R

This script requires that you first run Script_1.R, as it utilizes normalized counts data, an organized clinical data file, and a normalized intergenic TE matrix. The script performs GLMnet analysis, evaluating all different models depicted in Figure 2A. It also plots significant features from those models in Figure 2B. Additionally, it plots violin plots of normalized counts for four examples of significant features that differ between immune hot and cold, as shown in Figure 2C.

The script also performs a GLM test to assess the association between IKZF, Te score, and Immune types in the model, while adjusting for batch and histology. This analysis is reported in the manuscript. Furthermore, the script conducts conditional independence (mutual information) tests to identify causal relationships between TEs, IKZF1, and the immune-hot/-cold phenotype, as reported in the manuscript. Calculation of correlations between TE features (Supplemental Figure 6B)  & TE score calculation are also performed by this script.


## Script_3.R

This script requires that you first run Script_1.R and Script_2.R.

This script retrieves gene lists for immune pathways from published data and calculates TE-scores for significant TE families. It then performs  partial correlation and plotting of Figure 3A, and correlation and plotting between CD274, TE score and IKZF1 in figure Figure 3B,  and Kaplan-Meier plotting and analysis for Figure 3C (for IKZF1 and TE score).

## Script_1_TCGA.R

This script is the initial step in analyzing raw TCGA RNA sequencing data obtained from the REDISCOVERTE pipeline run. It normalizes the counts, filters out low counts, and performs immune deconvolution using the MCP counter. The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major clusters: immune hot and cold. The script also generates a heatmap of the two immune clusters and relevant clinical correlates, as shown in Figure 4A,B,C. Additionally, the script performs Kaplan-Meier analysis and plots for Immune hot/cold groups and Overall survival (OS), as seen in Figure 4D, and for IKZF1 high/low expression and OS survival for Figure 4F. Furthermore, the script conducts Cox Hazardous model analysis for overall survival and immune type for TCGA data.

## Script_2_TCGA.R

This script requires that you first run the Script_1_TCGA.R, as it utilizes normalized counts data, an organized Clinical TCGA file, and a normalized intergenic TE matrix from TCGA.

The script performs GLMnet analysis, evaluating all different models depicted in Supplemental Figure 8A. It also plots significant features from those models in Supplemental Figure 8B. Additionally, it plots violin plots of normalized counts for three examples of significant features that differ between immune hot and cold, as shown in Supplemental Figure 8C.


## Script_3_TCGA.R

This script requires that you first run the Script_1_TCGA.R and Script_2_TCGA.R

This script reads gene lists for immune pathways from published data (provided in data), calculates TE-scores for significant TE families, and performs partial correlation and plotting of Figure 9A, and correlation and plotting between CD274, TE score and IKZF1 in figure 9B. 
Additionally, it conducts Kaplan-Meier plotting and analysis for Figure 4E.


## AdditionalScripts
The AdditionalScripts directory includes scripts used for various analyses, such as:

Xcell_analysis.R: A script for processing RNA-seq counts and obtaining 64 cell types that are used for lymphoid and myeloid cell content calculation.

myeloid_and_lymphoid_content_calculation.R: Scripts for calculating myeloid and lymphoid content for each sample."

Supplemental_figure_4ABC.R ; Script for plotting purity for immune hot and cold, purity simulaiton, and person correlation between purity, lymphoid and myeloid cell content

Note: Purity estimates show in the figure, and disccused in the manuscript were retrieved by processing WES through TEMPO pipline (https://github.com/mskcc/tempo) as described in 
Materials and Methods section of the paper.

## HTML file reports
We also provide compiled HTML report files for each script, containing all the code and figures/tables presented in the manuscript. To access these reports, please download the HTM files for each script.

Study HTML reports: 

Script_1.html

Script_2.html

Script_3.html

Supplemental_figure_4ABC.html

TCGA HTML reports: 

Script_1_TCGA.html

Script_2_TCGA.html

Script_3_TCGA.html


## Contact
E-mail any questions to bradicm@mskcc.org


