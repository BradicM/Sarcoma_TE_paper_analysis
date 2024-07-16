## "Immune checkpoint inhibitor response in sarcomas associates with immune infiltrates and increased expression of transposable elements and viral response pathways" 

# Data description

The data directory includes the source data and code necessary to recreate all main and supplemental figures. The data directory has two subdirectories: 'study' and 'TCGA'.
The 'study' subdirectory features analysis of the clinical dataset collected at MSKCC, while the 'TCGA' subdirectory represents the data set and analysis that replicate our observations.

All necessary data files to reproduce all the figures are available in the 'data/study' and 'data/TCGA' folders

NOTE: study and TCGA folders do not include gene and Transposable element count matrices ("gene_expression counts" and "INTERGENIC_TE_normalized_expression_1052_repeats"
becaus of their size), and these files could be downloaded from Zenodo: https://zenodo.org/uploads/10313854 for TCGA data, and study data are avable with controlled access from 
dbgap acccession phs003284. 

# Code description 
All the costume code used to process RNA-seq data after intial processing by REDISCOVERTE can be found in the 'data/study/' folder, 'data/SAR038', 'data/TCGA/', 'data/study_on_trx' folder.




#  OUR STUDY ANALYSIS CODE


## Script_1.R 

This script performs a comprehensive analysis of our clinical raw RNA sequencing counts obtained from the REDISCOVERTE pipeline run. The script first normalizes and filters the counts, and then performs immune deconvolution using the MCP counter. The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major clusters: immune hot and cold.

The script generates several figures, including a heatmap of the immune clusters and relevant clinical correlates (Figure 1A,B,C), a Kaplan-Meier analysis and plot for Immune hot/cold groups and PFS (months) (Figure 1D), factor miner clustering of deconvoluted data (Supplemental Figure 1), a Cox Hazardous model analysis and plot (Supplemental Figure 2), and a matching of our immune hot and cold clusters with previously published Petitprez et.al clustering of immune cells (Supplemental Figure 3). Additionally, the script generates a heatmap of normalized intergenic TE counts with clinical correlates (Supplemental Figure 5).

The script also includes several supplementary figures, including a correlation between IKZF1 and B-cells (Supplemental Figure 6A) and number of reads for TE (Supplemental Figure 9A)
.

## Script_2.R

This script requires that you first run Script_1.R, as it utilizes normalized counts data, an organized clinical data file, and a normalized intergenic TE matrix. The script performs GLMnet analysis, evaluating all different models depicted in Figure 2A. It also plots significant features from those models in Figure 2B. Additionally, it plots violin plots of normalized counts for four examples of significant features that differ between immune hot and cold, as shown in Figure 2C.

The script also performs a GLM test to assess the association between IKZF, Te score, and Immune types in the model, while adjusting for batch and histology. This analysis is reported in the manuscript. Furthermore, the script conducts conditional independence (mutual information) tests to identify causal relationships between TEs, IKZF1, and the immune-hot/-cold phenotype, as reported in the manuscript. Calculation of correlations between TE features (Supplemental Figure 6B)  & TE score calculation are also performed by this script.


## Script_3.R

This script requires that you first run Script_1.R and Script_2.R.

This script retrieves gene lists for immune pathways from published data and calculates TE-scores for significant TE families. It then performs  partial correlation and plotting of Figure 3A, and correlation and plotting between CD274, TE score and IKZF1 in figure Figure 3B,  and Kaplan-Meier plotting and analysis for Figure 3C (for IKZF1 and TE score).




#   SARC038 COHORT ANALYSIS CODE


These are 3 scripts that are doing the same type of analysis that is described under OUR STUDY ANALYSIS CODE. The data set used here is SARC038 data that has been 
published before (#Reference paper and data https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1226445/full#supplementary-material
) and we use it here to confirm our observatinos. 
## Script_1.R 
## Script_2.R 
## Script_3.R 
These set of scripts plots Figure 5ABCDE, and Supplemental Figures 7, 8, 9A, 10. 
The order of running these scripts is the same as above, Script_1.R, Script_2.R,Script_3.R. 



#   TCGA ANALYSIS


## Script_1_TCGA.R

This script is the initial step in analyzing raw TCGA RNA sequencing data obtained from the REDISCOVERTE pipeline run. It normalizes the counts, filters out low counts, and performs immune deconvolution using the MCP counter. The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major clusters: immune hot and cold. The script also generates a heatmap of the two immune clusters and relevant clinical correlates, as shown in Figure 5A,B,C. Additionally, the script performs Kaplan-Meier analysis and plots for Immune hot/cold groups and Overall survival (OS), as seen in Figure 5D, and for IKZF1 high/low expression and OS survival for Figure 5F. Furthermore, the script conducts Cox Hazardous model analysis for overall survival and immune type for TCGA data.

## Script_2_TCGA.R

This script requires that you first run the Script_1_TCGA.R, as it utilizes normalized counts data, an organized Clinical TCGA file, and a normalized intergenic TE matrix from TCGA.

The script performs GLMnet analysis, evaluating all different models depicted in Supplemental Figure 12A. It also plots significant features from those models in Supplemental Figure 12B. Additionally, it plots violin plots of normalized counts for three examples of significant features that differ between immune hot and cold, as shown in Supplemental Figure 12C.


## Script_3_TCGA.R

This script requires that you first run the Script_1_TCGA.R and Script_2_TCGA.R

This script reads gene lists for immune pathways from published data (provided in data), calculates TE-scores for significant TE families, and performs partial correlation and plotting of Figure 13A, and correlation and plotting between CD274, TE score and IKZF1 in figure 13B. 
Additionally, it conducts Kaplan-Meier plotting and analysis for Figure 5 E.

This script also conducts Kaplan-Meier plotting and analysis for inindividual sarcima histology:  DDLP, UPS, and LMS  from TCGA data to determine TE score, IKZF1, and immune type association with PFS. 

Finally this script uses previously puvblished methoylation clustering data for DDLP, and correlates it with TE score. 


#  STUDY ON TREATMENT SAMPLE ANALYSIS


This script performs analysis of immune cell proportion clustering on baseline and on-treatment (immune therapy) data and determines which samples have switched/did not switch immune type and groups samples based on 4 categories:

1. Cold to hot switch
2. Hot to cold switch
3. Cold stable
4. Hot stable
   
Changes in IKZF and TE scores are then shown for each of these categories in baseline and on-treatment (Figure 6). Table 1 and Supplemental Figure 17 also show changes in immune cell proportions and associated statistics between baseline and on-treatment.

##The script that needs to be run for this analysis is Script_baseline_ontrx.R. 


## AdditionalScripts
The AdditionalScripts directory includes scripts used for various analyses, such as:

Xcell_analysis.R: A script for processing RNA-seq counts and obtaining 64 cell types that are used for lymphoid and myeloid cell content calculation.

myeloid_and_lymphoid_content_calculation.R: Scripts for calculating myeloid and lymphoid content for each sample."

Supplemental_figure_4ABC.R ; Script for plotting purity for immune hot and cold, purity simulaiton, and person correlation between purity, lymphoid and myeloid cell content

Note: Purity estimates show in the figure, and disccused in the manuscript were retrieved by processing WES through TEMPO pipline (https://github.com/mskcc/tempo) as described in 
Materials and Methods section of the paper.

## HTML file reports
We also provide compiled HTML report files for each script, containing all the code and figures/tables presented in the manuscript. To access these reports, please download the HTM files for each script.

##Study HTML reports: 

Script_1.html

Script_2.html

Script_3.html

Supplemental_figure_4ABC.html

##SARC038 html reports: 
Script_1.html

Script_2.html

Script_3.html


##TCGA HTML reports: 

Script_1_TCGA.html

Script_2_TCGA.html

Script_3_TCGA.html


## Contact
E-mail any questions to bradicm@mskcc.org or mb3188@gmail.com


