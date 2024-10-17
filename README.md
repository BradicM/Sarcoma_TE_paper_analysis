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


# System Requirements
Analysis requires only a standard computer with enough RAM to support the in-memory operations, and installation of already available R packages. 
Analysis wa performed on macOS Ventura 13.6.9, R version 4.3.1 (2023-06-16), and following R pacages versions were used: 

attached base packages:
splines,grid,stats4,stats,graphics,grDevices

other attached packages:
jtools_2.2.2,bnlearn_4.9.1,corrplot_0.92,ppcor_1.1,MASS_7.3-60.0.1,
dotwhisker_0.7.4,glmnet_4.1-8,Matrix_1.6-5,ggalt_0.4.0,gg.layers_0.1.1,
ggsignif_0.6.4,reshape2_1.4.4,splineTimeR_1.28.0,FIs_1.28.0,GeneNet_1.2.16,
fdrtool_1.2.17,longitudinal_1.1.13,corpcor_1.6.10,gtools_3.9.5,GSEABase_1.62.0,
graph_1.78.0,annotate_1.78.0,XML_3.99-0.16.1,AnnotationDbi_1.62.2,igraph_2.0.2,
tibble_3.2.1,GSVA_1.48.3,survminer_0.4.9,ggpubr_0.6.0,ggsci_3.0.0,
factoextra_1.0.7,FactoMineR_2.9,immunedeconv_2.1.0,EPIC_1.1.7,survival_3.5-7,
edgeR_3.42.4,limma_3.56.2,ggplot2_3.4.4,ComplexHeatmap_2.16.0,dplyr_1.1.4,
SummarizedExperiment_1.30.2,Biobase_2.60.0,GenomicRanges_1.52.0,GenomeInfoDb_1.36.4,IRanges_2.36.0,
S4Vectors_0.40.2,BiocGenerics_0.48.1,MatrixGenerics_1.14.0,matrixStats_1.2.0,tidyr_1.3.1,
data.table_1.15.0

loaded via a namespace (and not attached):
ggtext_0.1.2,fs_1.6.3,bitops_1.0-7,lubridate_1.9.3,insight_0.19.7,
ash_1.0-15,httr_1.4.7,RColorBrewer_1.1-3,doParallel_1.0.17,tools_4.3.1,
backports_1.4.1,utf8_1.2.4,R6_2.5.1,DT_0.31,HDF5Array_1.28.1,
mMCPcounter_1.1.0,mgcv_1.9-1,rhdf5filters_1.12.1,GetoptLong_1.0.5,withr_3.0.0,
prettyunits_1.2.0,gridExtra_2.3,preprocessCore_1.62.1,cli_3.6.2,exactRankTests_0.8-35,
Cairo_1.6-2,flashClust_1.01-2,sandwich_3.1-0,labeling_0.4.3,ComICS_1.0.4,
sass_0.4.8,mvtnorm_1.2-4,survMisc_0.5.6,readr_2.1.5,genefilter_1.82.1,
yulab.utils_0.1.4,commonmark_1.9.1,maps_3.4.2,readxl_1.4.3,rstudioapi_0.15.0,
RSQLite_2.3.5,gridGraphics_0.5-1,generics_0.1.3,shape_1.4.6,testit_0.13,
car_3.1-2,leaps_3.1,fansi_1.0.6,abind_1.4-5,terra_1.7-65,
lifecycle_1.0.4,scatterplot3d_0.3-44,multcomp_1.4-25,yaml_2.3.8,carData_3.0-5,
MCPcounter_1.2.0,ggstance_0.3.6,rhdf5_2.44.0,SparseArray_1.2.2,BiocFileCache_2.8.0,
blob_1.2.4,crayon_1.5.2,lattice_0.22-5,beachmat_2.18.0,cowplot_1.1.3,
KEGGREST_1.40.1,magick_2.8.2,pillar_1.9.0,knitr_1.45,boot_1.3-28.1,
rjson_0.2.21,estimability_1.4.1,codetools_0.2-19,gggrid_0.2-0,glue_1.7.0,
vctrs_0.6.5,png_0.1-8,cellranger_1.1.0,gtable_0.3.4,datawizard_0.9.1,
ggpp_0.5.6,cachem_1.0.8,ggpattern_1.0.1,xfun_0.42,S4Arrays_1.2.0,
coda_0.19-4,SingleCellExperiment_1.22.0,iterators_1.0.14,KMsurv_0.1-5,ellipsis_0.3.2,
TH.data_1.1-2,nlme_3.1-164,bit64_4.0.5,progress_1.2.3,filelock_1.0.3,
data.tree_1.1.0,bslib_0.6.1,maxstat_0.7-25,irlba_2.3.5.1,KernSmooth_2.23-22,
colorspace_2.1-0,DBI_1.2.1,tidyselect_1.2.0,emmeans_1.10.0,extrafontdb_1.0,
bit_4.0.5,compiler_4.3.1,curl_5.2.0,xml2_1.3.6,DelayedArray_0.28.0,
bayestestR_0.13.1,scales_1.3.0,proj4_1.0-14,multcompView_0.1-9,rappdirs_0.3.3,
stringr_1.5.1,digest_0.6.34,fftwtools_0.9-11,rmarkdown_2.25,XVector_0.42.0,
htmltools_0.5.7,pkgconfig_2.0.3,extrafont_0.19,sparseMatrixStats_1.14.0,highr_0.10,
dbplyr_2.4.0,fastmap_1.1.1,rlang_1.1.3,GlobalOptions_0.1.2,htmlwidgets_1.6.4,
DelayedMatrixStats_1.24.0,ggh4x_0.2.8,farver_2.1.1,jquerylib_0.1.4,zoo_1.8-12,
jsonlite_1.8.8,BiocParallel_1.36.0,BiocSingular_1.18.0,RCurl_1.98-1.14,magrittr_2.0.3,
polynom_1.4-1,ggplotify_0.1.2,GenomeInfoDbData_1.2.10,parameters_0.21.3,Rhdf5lib_1.22.1,
munsell_0.5.0,Rcpp_1.0.12,stringi_1.8.3,zlibbioc_1.48.0,plyr_1.8.9,
rtrend_0.1.5,parallel_4.3.1,ggrepel_0.9.5,Biostrings_2.68.1,pander_0.6.5,
gridtext_0.1.5,hms_1.1.3,circlize_0.4.15,locfit_1.5-9.8,markdown_1.12,
biomaRt_2.56.1,ScaledMatrix_1.10.0,evaluate_0.23,tzdb_0.4.0,foreach_1.5.2,
Rttf2pt1_1.3.12,purrr_1.0.2,km.ci_0.5-6,clue_0.3-65,rsvd_1.0.5,
broom_1.0.5,xtable_1.8-4,rstatix_0.7.2,memoise_2.0.1,cluster_2.1.6,
timechange_0.3.0,sva_3.48.0

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


#  ON TREATMENT ANALYSIS FROM OUR STUDY


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


