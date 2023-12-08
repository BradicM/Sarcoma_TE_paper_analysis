#Author Martina Bradic
#This script is the initial step in analyzing raw TCGA RNA sequencing data obtained from the REDISCOVERTE pipeline run.
#It normalizes the counts, filters out low counts, and performs immune deconvolution using the MCP counter. 
#The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major clusters: 
#immune hot and cold. The script also generates a heatmap of the two immune clusters and relevant clinical correlates, 
#as shown in Figure 4A,B,C. Additionally, the script performs Kaplan-Meier analysis and plots for Immune hot/cold groups 
#and Overall survival (OS), as seen in Figure 4D, and for IKZF1 high/low expression and OS survival for Figure 4F.
#Furthermore, the script conducts Cox Hazardous model analysis for overall survival and immune type for TCGA data.

library(data.table)
library(tidyr)
library("SummarizedExperiment")
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(edgeR)
library(survival)
library(immunedeconv)
library(FactoMineR)
library (factoextra)
library(ggsci)
library (factoextra)
library(survminer)

`%nin%` = Negate(`%in%`)



setwd("/Users/mbradic/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/!GDriveMigratedData/Nacev_project/TRANSPOSONS/FINAL_PAPER_ANALYSIS/PAPER_CODE_REPO/GITHUB_FOR_UPLOAD_FINAL/data/TCGA/")

#Set to the directory where scripts and all data are for analysis for TCGA data
#setwd("name of directory")

`%nin%` = Negate(`%in%`)

###################################################  GET GENE MATRIX FROM REDISCOVERTE QUANTIFICATION ######################

#Note: These are normalzied CPM log2 values which are obtained by the scaling method in edgeR, 
## These counts still need to be filtered for low expression and non-variable genes, and normalized using voom method
### prior to differential expresion and linear model analysis 

############ GET gene expression

gene_expression_mtrx<-read.table("gene_expression_mtrx")
colnames(gene_expression_mtrx)<-gsub(".","-",fixed=TRUE,colnames(gene_expression_mtrx))

#get batches info
TCGA_SARC_COHORT_BATCHES<-fread("TCGA_SARC_COHORT_BATCHES")

#this is a clinical file from the main TCGA sarcoma paper table S1
TCGA_SARC_clinical <- fread("TCGA_SARC_clinical.txt")
TCGA_SARC_clinical$patient<-gsub("-01","",fixed=TRUE, TCGA_SARC_clinical$`TCGA barcode`)


TCGA_SARC_clinical_ordered<-data.frame(colnames(gene_expression_mtrx))
colnames(TCGA_SARC_clinical_ordered)<-"patient"
RNAseq_pheno<-left_join(TCGA_SARC_clinical_ordered,TCGA_SARC_clinical)

all.equal(colnames(gene_expression_mtrx) ,RNAseq_pheno$patient)


#### Filter gene expression matrix filterByExpr

keep.exprs <- filterByExpr(gene_expression_mtrx, group=RNAseq_pheno$`short histo`)
gene_expression_mtrx_filtered <- gene_expression_mtrx[which(keep.exprs==TRUE),]
dim(gene_expression_mtrx_filtered)

#    after filtering 
#26010   190

### normalize it with voom
normalized_counts<-voom(gene_expression_mtrx_filtered)
normalized_counts<-normalized_counts$E



# remove genes that have no variance in expression values
#rows are samples columns are genes
#does not remove anything
normalized_counts_no_low_variance <- normalized_counts[which(apply(normalized_counts,1,var) > 0),]


##Subset batch file to get only those samples that we will analyze

TCGA_SARC_COHORT_BATCHES$patient<-TCGA_SARC_COHORT_BATCHES$`TCGA barcode`

#get only those that are in our analysis
TCGA_SARC_COHORT_BATCHES<-TCGA_SARC_COHORT_BATCHES[TCGA_SARC_COHORT_BATCHES$patient %in% RNAseq_pheno$patient, ]


#there are some duplicates in TCGA data set, let's find them 
TCGA_SARC_COHORT_BATCHES$patient[duplicated(TCGA_SARC_COHORT_BATCHES$patient)]

#check what they are
TCGA_SARC_COHORT_BATCHES[TCGA_SARC_COHORT_BATCHES$patient %in% c("TCGA-FX-A2QS","TCGA-K1-A3PN","TCGA-K1-A3PO","TCGA-K1-A42X","TCGA-SI-A71O","TCGA-VT-A80J"),]
#they seem to be comming from recurrent, normal or metastatis, and their batch ID of duplicates is always the same, 
#so we will select only primary tumors for this

TCGA_SARC_COHORT_BATCHES<-TCGA_SARC_COHORT_BATCHES[TCGA_SARC_COHORT_BATCHES$sample_type_name %in% "Primary Tumor",]

RNAseq_pheno<-left_join(RNAseq_pheno,TCGA_SARC_COHORT_BATCHES, by="patient")

all.equal(colnames(normalized_counts), RNAseq_pheno$patient)

normalized_counts<-normalized_counts[,colnames(normalized_counts) %in% RNAseq_pheno$patient]

#check for BATCH 
prcomp_values<-prcomp(t(normalized_counts))
var_explained <- prcomp_values$sdev^2/sum(prcomp_values$sdev^2)

lable_Prcomp<- data.frame(rownames(prcomp_values$x))
colnames(lable_Prcomp)<-"patient"
#lable_Prcomp$label<-gsub(".","_",lable_Prcomp$label, fixed=TRUE)
lable_Prcomp<-left_join(lable_Prcomp,RNAseq_pheno)



#Color by batch ID
ggplot(data.frame(prcomp_values$x), aes(x=PC1,y=PC2, label=RNAseq_pheno$batch_id, color=RNAseq_pheno$batch_id)) +
  geom_point(size = 3) +
  #geom_label(aes(fill = lable_Prcomp$Category), colour = "white", fontface = "bold")+
  theme_bw(base_size=32) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
  theme(legend.position="bottom")



###################################################  PERFORM IMMUNE DECONVOLUTION ######################

###contring with mcp_counter method     
res_mcp_counter = deconvolute(normalized_counts, "mcp_counter")

res_mcp_counter_table_2<-as.data.frame(res_mcp_counter %>% gather(sample, fraction, -cell_type))
res_mcp_counter_table_2_matrix<-spread(res_mcp_counter_table_2, cell_type,fraction)

#res_mcp_counter_table<-as.data.frame(res_mcp_counter)
res_mcp_counter_table<-dplyr::left_join(res_mcp_counter_table_2_matrix,RNAseq_pheno,by=c("sample"="patient"))
rownames(res_mcp_counter_table)<-res_mcp_counter_table$sample
res_mcp_counter_table<-res_mcp_counter_table[,-1]
res_mcp_counter_table_to_plot<-res_mcp_counter_table[,1:11]

mat2<-data.frame(res_mcp_counter_table_to_plot)

####Scale the values 
mat_immuno2<-t(scale(mat2))

all.equal(colnames(mat_immuno2),RNAseq_pheno$patient)
all.equal(rownames(mat2),RNAseq_pheno$patient)



#Merging all the STLMS, and ULMS into LMS

RNAseq_pheno$`short histo`<-gsub("STLMS" ,"LMS",RNAseq_pheno$`short histo`)
RNAseq_pheno$`short histo`<-gsub("ULMS" ,"LMS",RNAseq_pheno$`short histo`)


mat_immuno2_no_batch<-removeBatchEffect(mat_immuno2, RNAseq_pheno$batch_id)


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
#As shown by the increase in between group sum of squares, I really think the k-means consolidation provides quite a bit of stability here - 
#I ran HC without a consolidation procedure (using FactoMineR and ConsensusClusterPlus) and there were quite a bit of negative silhouette widths in the latter especially - which I think the consolidation procedure is addressing


#Lets export the cluster assignments so that we can compare features in more detail
#write.csv(res.hc[["data.clust"]],"[specify name of output file]")

cluster_assignments<-res.hc[["data.clust"]]



# Dendrogram, if dendrogram is needed

#fviz_dend(res.hc, 
#          cex = 0.7,                     # Label size
#          palette = "jco",               # Color palette see ?ggpubr::ggpar
#          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#          rect_border = "jco",           # Rectangle color
#          labels_track_height = 0.8      # Augment the room for labels
#)

# Individuals facor map, Supplemental Figure 7
fviz_cluster(res.hc, geom = "point", main = "Factor") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=10),axis.title.y = element_text(size=10), 
        axis.text.x = element_text( size=10),axis.text.y = element_text(size=10)) +
  scale_color_manual(values = c("1" = "red", "2" = "#0000FF"))



#To display quantitative variables that describe the most each cluster, type this:
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

cluster_assignments$Immune_subtypes<-ifelse(cluster_assignments$clust %in% 1, "Immune cold","Immune hot")


cluster_assignments$Immune_subtypes<-as.factor(cluster_assignments$Immune_subtypes)
cluster_assignments$patient<-rownames(cluster_assignments)
Ordered_pheno<-RNAseq_pheno
Ordered_pheno<-left_join(Ordered_pheno,cluster_assignments[,c("patient","Immune_subtypes")], by=c("patient"))
Ordered_pheno$`Abbreviation for Figures`<-Ordered_pheno$`short histo`
Ordered_pheno$batch <-Ordered_pheno$batch_id

table(Ordered_pheno$Immune_subtypes,Ordered_pheno$`Abbreviation for Figures`)

############################################################################################################################
#             Plotting heatmaps Figure 4ABC 
############################################################################################################################

mat_for_heatmap<-data.frame(t(mat_immuno2_no_batch))
mat_for_heatmap$patient<-rownames(mat_for_heatmap)
mat_for_heatmap_ordered<-data.frame(Ordered_pheno$patient)
colnames(mat_for_heatmap_ordered)<-"patient"
mat_for_heatmap_ordered<-left_join(mat_for_heatmap_ordered, mat_for_heatmap)
rownames(mat_for_heatmap_ordered)<-mat_for_heatmap_ordered$patient
mat_for_heatmap_ordered<-mat_for_heatmap_ordered[,c(-1)]


ann_for_2_types <- data.frame(Ordered_pheno[,c("short histo")]) 
rownames(ann_for_2_types)<-Ordered_pheno$patient
colnames(ann_for_2_types) <- c("Subtype") 

colours <- list(
  'Subtype'=c("DDLPS"= "#DEABCA","LMS"= "brown4","MFS"="#FF1493","UPS"= "#82D992")
  
)


colAnn <- HeatmapAnnotation(df = ann_for_2_types,
                            which = "col",
                            col = colours,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1,"mm"))



mat_for_heatmap_ordered<-scale(mat_for_heatmap_ordered)   # Z-score scaling by column. COlumns must be cells , rows should be samples

#replace "." in name of the cell type with space 
colnames(mat_for_heatmap_ordered)<-gsub("."," ",fixed=TRUE,colnames(mat_for_heatmap_ordered))


ht4AB =Heatmap(t(mat_for_heatmap_ordered), name = "mcp_counter_Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
             top_annotation =colAnn,show_column_dend = TRUE,show_row_dend = FALSE,show_column_names = FALSE)


#get expression values for genes that are part of PDl1 pathway and add them to heatmap 
mat_PDL1_related_genes_matrix<-normalized_counts[rownames(normalized_counts) %in% c("CTLA4","CD274","LAG3"),]

all.equal(colnames(mat_PDL1_related_genes_matrix),Ordered_pheno$patient)

mat_PDL1_related_genes_matrix_no_batch<-removeBatchEffect(mat_PDL1_related_genes_matrix, Ordered_pheno$batch_id,Ordered_pheno$`short histo`)

ht4C =Heatmap(mat_PDL1_related_genes_matrix_no_batch, name = "Epi_express",clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
             show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE)

#COMBINE HEATMAPS vertically  
FIG4_ABC = ht4AB %v% ht4C 
draw(FIG4_ABC)

############################################################################################################################

#             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR Immune hot/cold, Figure 4D and for IKZF1, figure 4F

############################################################################################################################

Ordered_pheno$OS_status<-ifelse(Ordered_pheno$`OS status` %in% "Dead",1,0)
Ordered_pheno$OS_status<-as.numeric(Ordered_pheno$OS_status)

# calculate months
Ordered_pheno$OS_months <-round(Ordered_pheno$`OS days`/30.417, digit=0)
Ordered_pheno$OS_months <-as.numeric(Ordered_pheno$OS_months)

#Model for immune hot/cold, Figure 4D
sfit <- survfit(Surv(OS_months, OS_status)~Immune_subtypes , data=Ordered_pheno)

ggsurvplot(sfit,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("#0000FF","red"),
           legend.labs=c("Immune-cold","Immune-hot"),font.legend=14,
           font.x=14,font.y=14,font.tickslab=14) 

fit <- coxph(Surv(OS_months, OS_status)~Immune_subtypes, data=Ordered_pheno)
summary(fit)

median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$OS_months)
median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$OS_months)



#Model for IKZF1 high/low, figure 4F
For_KLM_IKZF1<-data.frame(Ordered_pheno,normalized_counts[rownames(normalized_counts) %in% "IKZF1",])
colnames(For_KLM_IKZF1)[69]<-"IKZF1"
#Let's derive a cutpoint for IKZF1 expression 

IKZF21.surv_rnaseq.cut <- surv_cutpoint(
  data=For_KLM_IKZF1,
  time = "OS_months",
  event = "OS_status",
  variables = "IKZF1"
)
summary(IKZF21.surv_rnaseq.cut)


For_KLM_IKZF1$IKZF_high_low <-ifelse (For_KLM_IKZF1$IKZF1  > summary(IKZF21.surv_rnaseq.cut)[1]$cutpoint, "high expression","low expression")

sfit_IKZF1 <- survfit(Surv(OS_months, OS_status)~IKZF_high_low, data=For_KLM_IKZF1)

ggsurvplot(sfit_IKZF1,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
            legend.labs=c("high IKZF1","low IKZF1"),font.legend=14,
            font.x=14,font.y=14,font.tickslab=14) 




fit_IKZF1<- coxph(Surv(OS_months, OS_status)~IKZF1 + batch_id, data=For_KLM_IKZF1)
summary(fit_IKZF1)

############################################################################################################################

#             GET Trasposable element counts from REDISCOVERTE, normalize and visualize as heatmap

############################################################################################################################

### get TEs expression only for samples that we are interested in and only for LTR, LINE, SINE

#Get only Intergenic TEs

INTERGENIC_TE_normalized_expression_1052_repeats<-read.table("INTERGENIC_TE_normalized_expression_1052_repeats")
colnames(INTERGENIC_TE_normalized_expression_1052_repeats)<-gsub(".","-",fixed=TRUE,colnames(INTERGENIC_TE_normalized_expression_1052_repeats))

TE_names_order<-data.frame(rownames(INTERGENIC_TE_normalized_expression_1052_repeats))
colnames(TE_names_order)<-"repName"
TE_family_annotation<-read.table("TE_family_annotation")

TE_names_order<-left_join(TE_names_order, TE_family_annotation)
TE_names_order$repClass<-gsub("LTR?", "LTR", fixed=TRUE, TE_names_order$repClass)
TE_names_order$repClass<-gsub("SINE?", "SINE", fixed=TRUE, TE_names_order$repClass)



#### Filter gene expressio matrix filterByExpr
#genes are rows, samples are columns, use expression per group as minimun 

keep.exprs <- filterByExpr(INTERGENIC_TE_normalized_expression_1052_repeats, group=Ordered_pheno$`short histo`)
INTERGENIC_TE_normalized_expression_1052_repeats_filtered <- INTERGENIC_TE_normalized_expression_1052_repeats[which(keep.exprs==TRUE),]
dim(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)
###There 1002 that are left here to be analysed

### normalize it with voom
normalized_intergenic_TE<-voom(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)
normalized_TE_counts<-normalized_intergenic_TE$E

###Continue analysis with normalized data
INTERGENIC_TE_normalized_expression_1052_repeats<-data.frame(t(normalized_TE_counts))

INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed<-t(INTERGENIC_TE_normalized_expression_1052_repeats %>% mutate_at(colnames(INTERGENIC_TE_normalized_expression_1052_repeats), scale))
colnames(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed)<-rownames(INTERGENIC_TE_normalized_expression_1052_repeats)

#remove those with NAs, these are 3 TEs acros all the samples only after we transfor with Z-score
INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed<-INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed[!rowSums(is.na(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed)) > 0,]

#checj is Phenotypes order match to Te matrix order so I can use phenotyes for heatmap 
all.equal(Ordered_pheno$patient,rownames(INTERGENIC_TE_normalized_expression_1052_repeats))


INTERGENIC_TE_normalized_expression_1052_repeats<-t(INTERGENIC_TE_normalized_expression_1052_repeats)


