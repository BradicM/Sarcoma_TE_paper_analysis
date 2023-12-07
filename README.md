## "Immune checkpoint inhibitor response in sarcomas associates with immune infiltrates and increased expression of transposable elements and viral response pathways" 

## Data and code
Data directory includes source data and code necessary to recreate all main and supplemental figures.
Data directory has two subdirectories; study, and TCGA featuring analysis of the clinical dataset collected at MSKCC,
and TCG represents data set and anysis that replicate our observations


Data files necessery to reproduce all the figures are avavible in data folder and include: 

Scripts to process RNA-seq are found in data/study/Scripts/ folder or data/TCGA/Scripts/ folder and include

```
Script_1.R

```

This is the initial script that reads in raw RNA seq counts obtained from REDISCOVERTE pipeline, normalizes counts, filters low counts, 
and performs immune deconvolution using MCP counter, and then FactoMineR clustering of the MCP counter deconvoluted immune cell proportions
to establish 2 major clusters; immune hot and cold. The first part of the script plots heatmap of those immune 
clusters and relevant clinical correlates; Figure 1A,B,C
This scripts also performs following analysis and plots Figures:
Kaplan-Meier analysis and plots for Immune hot/cold groups and PFS (months); Figure 1D
Cox Hazardous model analysis and plot Supplemental Figure 2 
Matching of our immune hot and cold clusters with previously published Petitprez et.al clustering of immune cells Supplemental Figure 3
Supplemental Figure 6A, correlation between IKZF1 and B-cells
The second part of the script reads in intergenic TE counts, normalizes and filters them, and plots those normalized data heatmap with clinical correlates; Supplemental figure 5

```
Script_2.R
```


```
Script_3.R
```

```
Script_1_TCGA.R
```

```
Script_2_TCGA.R
```

```
Script_3_TCGA.R
```

## AdditionalScripts
AdditionalScripts directory includes scripts used for for following analysis: 

```
Script for processing RNA-seq counts and to obtain 64 cell types that are used for lymphoid and myeloid cell content calculation
Xcell_analysis.R

Scripts for calcualtion of myeloid and lymphoid content for each sample
myeloid_and_lymphoid_content_calculation.R

```


## Contact
E-mail any questions to bradicm@mskcc.org


