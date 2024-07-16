#Author Martina Bradic
#This script performs a comparison between baseline and on-trx paires. 
#First we cluster all the samples baseline and on-trx and determine their immune type; Supplemental Figure 17
#we then compare changes in immune proportion,  TE expression and IKZF1 expression between pairs pre and post treatment 
# this script generates Figure 6, and Table 1, as well as all associated statistics


library(data.table)
library(tidyr)
library("SummarizedExperiment")
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(edgeR)
library(survival)
library(survminer)
library(ggsci)
library (factoextra)
library(FactoMineR)
library (factoextra)
library(immunedeconv)
library(GSVA)

`%nin%` = Negate(`%in%`)


###################################################  GET GENE MATRIX FROM REDISCOVERTE QUANTIFICATION ######################

#Note: These are normalzied CPM log2 values which are obtained by the scaling method in edgeR, 
## These counts still need to be filtered for low expression and non-variable genes, and normalized using voom method
### prior to differential expression and linear model analysis 

Clinical_phenotypes_ontrx<-read.table("Clinical_phenotypes_ontrx")
gene_expression_mtrx_ontrx<-read.table("gene_expression_mtrx_ontrx")
colnames(gene_expression_mtrx_ontrx)<-gsub(".","-",fixed=TRUE,colnames(gene_expression_mtrx_ontrx))

names_ontrx<-data.frame(colnames(gene_expression_mtrx_ontrx))
colnames(names_ontrx)<-"Pathology_Accession_nb"
names_ontrx<-left_join(names_ontrx, Clinical_phenotypes_ontrx)
all.equal(names_ontrx$Pathology_Accession_nb,colnames(gene_expression_mtrx_ontrx))

colnames(gene_expression_mtrx_ontrx)<-names_ontrx$CMO_ID
colnames(gene_expression_mtrx_ontrx)<-paste(colnames(gene_expression_mtrx_ontrx), sep="_","Ontrx")

gene_expression_mtrx_baseline<-read.table("gene_expression_mtrx_baseline")
colnames(gene_expression_mtrx_baseline)<-gsub(".","-",fixed=TRUE,colnames(gene_expression_mtrx_baseline))

colnames(gene_expression_mtrx_baseline)<-paste(colnames(gene_expression_mtrx_baseline), sep="_","base")


all.equal(rownames(gene_expression_mtrx_ontrx),rownames(gene_expression_mtrx_baseline))

gene_expression_mtrx<-cbind(gene_expression_mtrx_baseline,gene_expression_mtrx_ontrx)

#get clinical data
Clinical_phenotypes<-read.table("Clinical_phenotypes_combined.txt", sep="\t", header=TRUE)


#order clinical file based on gene expresion 
Ordered_pheno<-data.frame(colnames(gene_expression_mtrx))
colnames(Ordered_pheno)<-"CMO_ID"
Ordered_pheno<-left_join(Ordered_pheno,Clinical_phenotypes)
Ordered_pheno$`Abbreviation for Figures`<-Ordered_pheno$Abbreviation.for.Figures

#perform filtering
keep.exprs <- filterByExpr(gene_expression_mtrx, group=Ordered_pheno$`Abbreviation for Figures`)
gene_expression_mtrx_filtered <- gene_expression_mtrx[which(keep.exprs==TRUE),]
dim(gene_expression_mtrx_filtered)

#    after filtering 
#29432    113

### normalize it with voom
normalized_counts<-voom(gene_expression_mtrx_filtered)
normalized_counts<-normalized_counts$E


###contring with mcp_counter mehtod     
res_mcp_counter = deconvolute(normalized_counts, "mcp_counter")

res_mcp_counter_table_2<-as.data.frame(res_mcp_counter %>% gather(sample, fraction, -cell_type))
res_mcp_counter_table_2_matrix<-spread(res_mcp_counter_table_2, cell_type,fraction)


res_mcp_counter_table<-dplyr::left_join(res_mcp_counter_table_2_matrix,Ordered_pheno,by=c("sample"="CMO_ID"))
rownames(res_mcp_counter_table)<-res_mcp_counter_table$sample
res_mcp_counter_table<-res_mcp_counter_table[,-1]
res_mcp_counter_table_to_plot<-res_mcp_counter_table[,1:11]

mat2<-data.frame(res_mcp_counter_table_to_plot)

####Scale the values 
mat_immuno2<-t(scale(mat2))


all.equal(colnames(mat_immuno2),Ordered_pheno$CMO_ID)

#need to reorder
correct_order_for_batch<-data.frame(colnames(mat_immuno2))
colnames(correct_order_for_batch)<-"CMO_ID"
correct_order_for_batch<-left_join(correct_order_for_batch,Ordered_pheno)
all.equal(colnames(mat_immuno2),correct_order_for_batch$CMO_ID)

mat_immuno2_no_batch<-removeBatchEffect(mat_immuno2, correct_order_for_batch$batch,correct_order_for_batch$`Abbreviation for Figures`)


##################################################### DETERMINE CLUSTERS FROM MCP counter using FactoMineR analysis ###############
#set seed for reproducibilty
set.seed(2021)

#We will apply agglomerative HC using Ward's method/Euclidean distance across a range of k from 2-12; we will have FactoMineR suggest the optimal k which is does based on within-cluster sum of squares (i.e. inertia);  the initial partition suggested by HC is followed by a k-means procedure to consolidate cluster membership
res.hc <- HCPC(as.data.frame(t(mat_immuno2_no_batch)), 
               nb.clust = -1,      #Allow FactoMineR to suggest optimal partition (minimize within cluster sum of sq)
               consol=TRUE,        #Perform k-means consolidation after initial HC partition
               iter.max = 10,      #Iterations for the k-means consolidation
               min = 2,            #Min k
               max = 12,           #Max k
               metric="euclidean", #Distance metric
               method="ward",      #Ward's method
               graph = FALSE)

#Lets evaluate between group sum of squares (inertia) before and after k-means consolidation - it should go up as we are maximizing separation
res.hc$call$bw.before.consol
res.hc$call$bw.after.consol

#Now that we've performed the clustering, lets visualize the HC dendrogram using factoextra, will use the lancet color palette from ggsci package 
#Note that this reflects the HC partition prior to the k-means consolidation


#Lets export the cluster assignments so that we can compare features in more detail
#write.csv(res.hc[["data.clust"]],"[specify name of output file]")

cluster_assignments<-res.hc[["data.clust"]]


# If you want to plot dendrogram 
#fviz_dend(res.hc, 
#          cex = 0.7,                     # Label size
#          palette = "jco",               # Color palette see ?ggpubr::ggpar
#          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#          rect_border = "jco",           # Rectangle color
#          labels_track_height = 0.8      # Augment the room for labels
#)

# Individuals facor map, This is Supplemental Figure 1
fviz_cluster(res.hc, geom = "point", main = "Factor") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10), 
        axis.text.x = element_text( size=10),axis.text.y = element_text(size=10)) +
  scale_color_manual(values = c("1" = "red", "2" = "#0000FF"))



#Supplemental table 2: Display quantitative variables that describe the most each cluster
res.hc$desc.var$quanti

#show principal dimensions that are the most associated with clusters, type this:
res.hc$desc.axes$quanti

#Finally, representative individuals of each cluster can be extracted as follow:
res.hc$desc.ind$para

Cluster_1_individuals_contrib<-unlist(res.hc$desc.ind$para[1])
Cluster_1_individuals_contrib<-names(Cluster_1_individuals_contrib)
Cluster_1_individuals_contrib<- gsub("1.","", Cluster_1_individuals_contrib)
Cluster_2_individuals_contrib<-unlist(res.hc$desc.ind$para[2])
Cluster_2_individuals_contrib<-names(Cluster_2_individuals_contrib)
Cluster_2_individuals_contrib<- gsub("2.","", Cluster_2_individuals_contrib)

#match orders of samples
all.equal(rownames(cluster_assignments), colnames(mat_immuno2))
all.equal(rownames(cluster_assignments), colnames(t(mat2)))
all.equal(rownames(cluster_assignments), colnames(mat_immuno2_no_batch))

#call cluster names
cluster_assignments$Immune_subtypes<-ifelse(cluster_assignments$clust %in% 1, "Immune cold","Immune hot")
cluster_assignments$CMO_ID<-rownames(cluster_assignments)


#Assign cluster names into phenotyo phile
Ordered_pheno<-left_join(Ordered_pheno,cluster_assignments[,c("CMO_ID","Immune_subtypes")])





############################################################################################################################

#             GET Trasposable element counts from REDISCOVERTE, normalize 
############################################################################################################################


##################################.   GET Transportable element counts from REDISCOVER TE pipeline and prepare for analysis  ##########################################
### get TEs expression only for LTR, LINE, SINE elements

#GEt only Intergenic baseline
INTERGENIC_TE_normalized_expression_1052_repeats_baseline<-read.table("INTERGENIC_TE_normalized_expression_1052_repeats_baseline")
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_baseline)<-gsub(".","-",fixed=TRUE,colnames(INTERGENIC_TE_normalized_expression_1052_repeats_baseline))


#read intergenic TEs ontrx

INTERGENIC_TE_normalized_expression_1052_repeats_ontrx<-read.table("INTERGENIC_TE_normalized_expression_1052_repeats_ontrx")
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx)<-gsub(".","-",fixed=TRUE,colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx))

dim(INTERGENIC_TE_normalized_expression_1052_repeats_baseline)
dim(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx)

#all 1052 are here 
sharedTEs<-intersect(rownames(INTERGENIC_TE_normalized_expression_1052_repeats_baseline),rownames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx))


#MERGE THEM ALL TOGETHER 

#renames sample names for baseline and ontrx first
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_baseline)<-paste(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_baseline), sep="_","base")


names_ontrx_TE<-data.frame(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx))
colnames(names_ontrx_TE)<-"Pathology_Accession_nb"
names_ontrx_TE<-left_join(names_ontrx_TE, Clinical_phenotypes_ontrx)
all.equal(names_ontrx_TE$Pathology_Accession_nb,colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx))

colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx)<-names_ontrx_TE$CMO_ID
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx)<-paste(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_ontrx), sep="_","Ontrx")



INTERGENIC_TE_normalized_expression_1052_repeats<-cbind(INTERGENIC_TE_normalized_expression_1052_repeats_baseline,INTERGENIC_TE_normalized_expression_1052_repeats_ontrx)



#Match TEs with their familes 
TE_names_order<-data.frame(rownames(INTERGENIC_TE_normalized_expression_1052_repeats))
colnames(TE_names_order)<-"repName"
TE_family_annotation<-read.table("TE_family_annotation")
TE_names_order<-left_join(TE_names_order, TE_family_annotation)
TE_names_order$repClass<-gsub("LTR?", "LTR", fixed=TRUE, TE_names_order$repClass)
TE_names_order$repClass<-gsub("SINE?", "SINE", fixed=TRUE, TE_names_order$repClass)


#### Filter gene TE matrix filterByExpr
#TEs are rows, samples are columns, use expression per group as a minimum 

#check the order of samples
all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats), Ordered_pheno$CMO_ID)

keep.exprs_TE <- filterByExpr(INTERGENIC_TE_normalized_expression_1052_repeats, group=Ordered_pheno$`Abbreviation for Figures`)
INTERGENIC_TE_normalized_expression_1052_repeats_filtered <- INTERGENIC_TE_normalized_expression_1052_repeats[which(keep.exprs_TE==TRUE),]
dim(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)

all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_filtered), Ordered_pheno$CMO_ID)

###There are 844 TEs that are left here to be analysed

### normalize it with voom
normalized_intergenic_TE<-voom(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)
normalized_TE_counts<-normalized_intergenic_TE$E

###Continue analysis with TE normalized data
INTERGENIC_TE_normalized_expression_1052_repeats<-data.frame(t(normalized_TE_counts))

INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed<-t(INTERGENIC_TE_normalized_expression_1052_repeats %>% mutate_at(colnames(INTERGENIC_TE_normalized_expression_1052_repeats), scale))
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed)<-rownames(INTERGENIC_TE_normalized_expression_1052_repeats)

#remove those with NAs, these are 3 TEs acros all the samples only after we transfor with Z-score
INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed<-INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed[!rowSums(is.na(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed)) > 0,]

INTERGENIC_TE_normalized_expression_1052_repeats<-t(INTERGENIC_TE_normalized_expression_1052_repeats)

#check if Prototypes order match to TE matrix order
all.equal(Ordered_pheno$CMO_ID,colnames(INTERGENIC_TE_normalized_expression_1052_repeats))

##################### PLOT INTERGENIC NORMALIZED EXPRESSION OVER TWO CLSUTERS

INTERGENIC_TE_normalized_expression_1052_repeats_no_batch<-removeBatchEffect(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed, Ordered_pheno$batch,Ordered_pheno$`Abbreviation for Figures`)

#TE_intergenic_for_heatmap<-t(scale(INTERGENIC_TE_normalized_expression_1052_repeats))

TE_intergenic_for_heatmap<-INTERGENIC_TE_normalized_expression_1052_repeats_no_batch



#For_heatmap_TE_annotation<-readRDS("~/juno/work/ccs/badicm/Project_NACEV/RNAseq/REdiscoverTE/REdiscoverTE/ANALYSIS_FINAL/Rollup_analysis/ALL_INCLUDING_DICSON_projects/RESULTS/RE_intergenic_2_counts_normalized.RDS")
TE_intergenic_for_heatmap_row_anno<-data.frame(rownames(TE_intergenic_for_heatmap))
colnames(TE_intergenic_for_heatmap_row_anno)<-"repName"

# Fix some TE family names, replace "." with "-" 
TE_intergenic_for_heatmap_row_anno$repName<-gsub("CR1.","CR1-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub(".int","-int",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("ERV3.","ERV3-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("ERVL.","ERVL-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("hAT.","hAT-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("HERV.","HERV-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("HUERS.","HUERS-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("L1PA15.","L1PA15-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("L2.","L2-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("MamGypsy2.","MamGypsy2-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)
TE_intergenic_for_heatmap_row_anno$repName<-gsub("ORSL.","ORSL-",TE_intergenic_for_heatmap_row_anno$repName,fixed = TRUE)


TE_intergenic_for_heatmap_row_anno<-left_join(TE_intergenic_for_heatmap_row_anno,TE_family_annotation)
TE_intergenic_for_heatmap_row_anno$repClass<-gsub("LTR?","LTR",fixed = TRUE,TE_intergenic_for_heatmap_row_anno$repClass)
TE_intergenic_for_heatmap_row_anno$repClass<-gsub("SINE?","SINE",fixed = TRUE,TE_intergenic_for_heatmap_row_anno$repClass)

#recheck if it was replaced
TE_intergenic_for_heatmap_row_anno[is.na(TE_intergenic_for_heatmap_row_anno$repClas),]


###Create annotations and plot Supplemental Figure 5

#Annotation for columns
ann_for_2_types <- data.frame(Ordered_pheno[,c("Response","Abbreviation for Figures", "Sample_Timepoint")]) 
rownames(ann_for_2_types)<-Ordered_pheno$"CMO_ID"
colnames(ann_for_2_types) <- c("Response","Subtype","Sample_Timepoint") 



colours_responders <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  'Subtype'=c ("ANGS"="#DD67C0" ,"CHS"="#A848E0","DDLS"= "#DEABCA","LMS"= "#F4A460","OS"= "#8881D2" ,"SARCNOS"="#8DB6D1","UPS"= "#82D992","EHE"="#FFFF00","LPS"="#0000FF","MFS"="#B22222","Other"="#FF0000","SBRC"="#85E0D6","ASPS"="grey"),
  
  'Sample_Timepoint'=c("Baseline"="blue",'Ontrx'="red")
)


colAnn2 <- HeatmapAnnotation(df = ann_for_2_types,
                             which = "col",
                             col = colours_responders,
                             annotation_width = unit(c(1, 4), "cm"),
                             gap = unit(1,"mm"))


#Annotation for rows

colours_rows <- list(
  'repClass' = c("DNA"="orange","SINE"="darkpurple" , "LINE"= "darkblue","LTR"="darkgreen")
)

rowAnn <- HeatmapAnnotation(df = TE_intergenic_for_heatmap_row_anno$repClass,
                            which = "row",
                            col = colours_rows,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1,"mm"))




################################.#######################################   CALCALATE TE scores  ########################################


TE_matrix<-INTERGENIC_TE_normalized_expression_1052_repeats
TE_matrix<-t(TE_matrix)
#Check if the ORder pheno is the same sample order with TE_matrix
all.equal(rownames(TE_matrix) , Ordered_pheno$CMO_ID)



INTERGENIC_TE_normalized_expression_1052_repeats[rownames(INTERGENIC_TE_normalized_expression_1052_repeats) %in% c("MER57F","Tigger17a","MER61F","MER45A","LTR104_Mam","HERVL74.int"),]


#calculate TE score for significant Te features 
TE_features<-c("Tigger17a","MER61F","MER45A","LTR104_Mam","HERVL74.int")

#geneSets<-list()
#geneSets$TE_features<-as.character(TE_features)
#geneSets<-lapply(geneSets, function(z){ z[!is.na(z) & z != ""]})

TE_genes_matrix<-as.matrix(TE_matrix[,colnames(TE_matrix) %in% TE_features])
TE_genes_matrix_df<-data.frame(TE_genes_matrix)
TE_genes_matrix_df$CMO_ID<-rownames(TE_genes_matrix_df)

TE_genes_matrix_df<-left_join(TE_genes_matrix_df, Ordered_pheno)


TE_features<-data.frame(TE_features)

#calculate TE scores
geneSets<-as.list(TE_features)
geneSets$TE_features<-as.character(geneSets$TE_features)
geneSets<-lapply(geneSets, function(z){ z[!is.na(z) & z != ""]})

TE_genes_matrix<-as.matrix(TE_matrix[,colnames(TE_matrix) %in% as.character(TE_features$TE_features)])

gsva_for_TE__genes <-gsva(t(TE_genes_matrix),geneSets,method="zscore", kcdf="Gaussian",verbose=FALSE)

#check for same order
all.equal(colnames(gsva_for_TE__genes),Ordered_pheno$CMO_ID)

TE_genes_matrix_df<- data.frame(t(gsva_for_TE__genes), Ordered_pheno)



##############################################     Function to plot regression.  ############################################## 
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


################################################.  COMPARISON BETWEEN PAIRS for immune cells, TEs, TLS and IKZF1 #########################################################
#load function to run this analysis

source("Analyze_individual_groups_function.R")


###Run only one comparison at the the time, and continue to "RUN PAIRS COMPARISON" section and run till the end

##COMPARISON 1
#get pairs that go from cold to hot
Cold_to_hot_switch<-fread("Cold_to_hot_switch.txt",sep="\t", header=TRUE)
get_individual_pairs_comparison_calculations(Cold_to_hot_switch)

#######################

##COMPARISON 2
#get pairs that go from hot to cold

Hot_to_cold_switch<-fread("Hot_to_cold_switch.txt",sep="\t", header=TRUE)
get_individual_pairs_comparison_calculations(Hot_to_cold_switch)
#######################


##COMPARISON 3
#get pairs that go from hot to cold
NO_switch_samples<-fread("NO_switch_samples.txt",sep="\t", header=TRUE)
NO_switch_samples_cold<-NO_switch_samples[NO_switch_samples$`After treatment` %in% "Immune cold",]
get_individual_pairs_comparison_calculations(NO_switch_samples_cold)
#######################


##COMPARISON 4
NO_switch_samples_hot<-NO_switch_samples[NO_switch_samples$`After treatment` %in% "Immune hot",]
get_individual_pairs_comparison_calculations(NO_switch_samples_hot)


################################################.  COMPARISON BETWEEN PAIRS for 24 immune pathways #########################################################
source("Analyze_individual_groups_for_24_pathways.R")

Gene_lists_for_ddGSEA<-fread("Gene_lists_for_ddGSEA.txt",header=TRUE)

C_GAS_gene_list_KEGG<-fread("C_GAS_gene_list_KEGG")
C_GAS_gene_list_KEGG<-C_GAS_gene_list_KEGG[-c(1),]
C_GAS_gene_list_KEGG<-data.frame(C_GAS_gene_list_KEGG)
colnames(C_GAS_gene_list_KEGG)<-"genes"



Gene_lists_for_ddGSEA<-Gene_lists_for_ddGSEA %>% 
  mutate(Genes=strsplit(Genes, ",")) %>% 
  unnest(Genes)

Gene_list<-split.data.frame(Gene_lists_for_ddGSEA[-1],Gene_lists_for_ddGSEA$`Gene Signature`)
geneSets<-lapply(Gene_list,"[[","Genes")

#add c-gas
geneSets$C_GAS_gene<-C_GAS_gene_list_KEGG$genes

##NOW calculate ssGSES for each sample using normalized counts 
gsva_24_gene_lists<-gsva(as.matrix(normalized_counts),geneSets,method="gsva",kcdf="Gaussian" ,verbose=FALSE)

all.equal(colnames(gsva_24_gene_lists),Ordered_pheno$CMO_ID)


gsva_24_gene_lists_df<-data.frame(t(gsva_24_gene_lists))
gsva_24_gene_lists_df$CMO_ID<-rownames(gsva_24_gene_lists_df)


get_individual_pairs_comparison_calculations_24_pathways(Cold_to_hot_switch)
get_individual_pairs_comparison_calculations_24_pathways(Hot_to_cold_switch)
get_individual_pairs_comparison_calculations_24_pathways(NO_switch_samples_cold)
get_individual_pairs_comparison_calculations_24_pathways(NO_switch_samples_hot)


#rmarkdown::render("Script_baseline_ontrx.R")
