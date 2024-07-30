#Author Martina Bradic
#This script performs a comprehensive analysis of our clinical raw RNA sequencing counts obtained from the REDISCOVERTE
#pipeline run. The script analyses SARC038 validation cohort
#Reference paper and data https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1226445/full#supplementary-material
###RNA sequencing libraries were constructed using the KAPA RNA Hyper kit to generate a total RNA library, which was further captured using the Twist Core Exome probe set which targets primarley exons, thus 
# not a high number of TE related genes are expected. 

#Scripts creates following figures
#Figure 4 ABCD
#Supplementary Figure 7 
#Supplementary Figure 9B

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
library(ggfortify)

`%nin%` = Negate(`%in%`)



###################################################  GET GENE MATRIX FROM REDISCOVERTE QUANTIFICATION ######################


#Note: These are normalzied CPM log2 values which are obtained by the scaling method in edgeR, 
## These counts still need to be filtered for low expression and non-variable genes, and normalized using voom method
### prior to differential expression and linear model analysis 

############ GET gene expression
gene_expression_mtrx<-read.table("gene_expression_mtrx")
rownames(gene_expression_mtrx)<-gene_expression_mtrx[,1]
gene_expression_mtrx<-gene_expression_mtrx[,-c(1)]

colnames(gene_expression_mtrx)<-gsub(".","-",fixed=TRUE,colnames(gene_expression_mtrx))

#get clinical data
Clinical_phenotypes<-fread("Clinical_phenotypes.txt")


#order clinical file based on gene expresion 
Ordered_pheno<-data.frame(colnames(gene_expression_mtrx))
colnames(Ordered_pheno)<-"Sample_ID"
Ordered_pheno<-left_join(Ordered_pheno,Clinical_phenotypes)
dim(Ordered_pheno)
#remove samples for which we have missing data

Ordered_pheno<-Ordered_pheno[! is.na(Ordered_pheno$NOI_ID),]
Ordered_pheno<-Ordered_pheno[Ordered_pheno$Timepoint %in% "baseline",]

#just take baseline samples
Ordered_pheno
dim(Ordered_pheno)
gene_expression_mtrx<-gene_expression_mtrx[,colnames(gene_expression_mtrx) %in% Ordered_pheno$Sample_ID,]
dim(gene_expression_mtrx)
#total of 36 samples
Ordered_pheno$`Abbreviation for Figures`<-Ordered_pheno$Diagnosis
Ordered_pheno$Response<-Ordered_pheno$`Best response`
Ordered_pheno$Response<-gsub("CR","PR",Ordered_pheno$Response)
Ordered_pheno$Response<-gsub("PR","CR/PR",Ordered_pheno$Response)

#perform filtering
keep.exprs <- filterByExpr(gene_expression_mtrx, group=Ordered_pheno$`Abbreviation for Figures`)
gene_expression_mtrx_filtered <- gene_expression_mtrx[which(keep.exprs==TRUE),]
dim(gene_expression_mtrx_filtered)

#    after filtering 
#23871    36

### normalize it with voom
normalized_counts<-voom(gene_expression_mtrx_filtered)
normalized_counts<-normalized_counts$E


#create a PCA plot
pca_res<-prcomp(t(normalized_counts))

all.equal(Ordered_pheno$Sample_ID,rownames(pca_res$x))
rownames(pca_res$x)<-Ordered_pheno$NOI_ID

autoplot(pca_res, data = Ordered_pheno, colour = 'Diagnosis',label = TRUE, label.size = 3)



###contring with mcp_counter mehtod     
res_mcp_counter = deconvolute(normalized_counts, "mcp_counter")

res_mcp_counter_table_2<-as.data.frame(res_mcp_counter %>% gather(sample, fraction, -cell_type))
res_mcp_counter_table_2_matrix<-spread(res_mcp_counter_table_2, cell_type,fraction)


res_mcp_counter_table<-dplyr::left_join(res_mcp_counter_table_2_matrix,Ordered_pheno,by=c("sample"="Sample_ID"))
rownames(res_mcp_counter_table)<-res_mcp_counter_table$sample
res_mcp_counter_table<-res_mcp_counter_table[,-1]
res_mcp_counter_table_to_plot<-res_mcp_counter_table[,1:11]

mat2<-data.frame(res_mcp_counter_table_to_plot)

####Scale the values 
mat_immuno2<-t(scale(mat2))


all.equal(colnames(mat_immuno2),Ordered_pheno$Sample_ID)

#need to reorder
correct_order_for_batch<-data.frame(colnames(mat_immuno2))
colnames(correct_order_for_batch)<-"Sample_ID"
correct_order_for_batch<-left_join(correct_order_for_batch,Ordered_pheno)
all.equal(colnames(mat_immuno2),correct_order_for_batch$Sample_ID)


#correct for histology

mat_immuno2_no_batch<-removeBatchEffect(mat_immuno2,correct_order_for_batch$`Abbreviation for Figures`)


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

# Individuals facor map, Supplementary Figure S7
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
cluster_assignments$Sample_ID<-rownames(cluster_assignments)


#Assign cluster names into phenotyo phile
Ordered_pheno<-left_join(Ordered_pheno,cluster_assignments[,c("Sample_ID","Immune_subtypes")])



########################  PLOTTING HEATMAPS ; Figure 4 ##########################################


#we remove batch before
mat_for_heatmap<-data.frame(t(mat_immuno2_no_batch))

mat_for_heatmap$Sample_ID<-rownames(mat_for_heatmap)
mat_for_heatmap_ordered<-data.frame(Ordered_pheno$Sample_ID)
colnames(mat_for_heatmap_ordered)<-"Sample_ID"
mat_for_heatmap_ordered<-left_join(mat_for_heatmap_ordered, mat_for_heatmap)
rownames(mat_for_heatmap_ordered)<-mat_for_heatmap_ordered$Sample_ID
mat_for_heatmap_ordered<-mat_for_heatmap_ordered[,c(-1)]

ann_for_2_types <- data.frame(Ordered_pheno[,c("Response","Abbreviation for Figures")]) 
rownames(ann_for_2_types)<-Ordered_pheno$Sample_ID
colnames(ann_for_2_types) <- c("Response","Subtype") 

colours <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  
  'Subtype'=c ( "DDLPS"= "#DEABCA","UPS"= "#82D992", "LS"= "#F4A460","OS"= "#8881D2" , "CS"="green","ES"="purple","SS"="magenta")
                
#synovial sarcomas (SS), four leiomyosarcomas (LMS), two dedifferentiated liposarcomas (DDLPS), three undifferentiated pleomorphic sarcomas (UPS); and 18 patients with bone sarcomas, including six Ewing sarcomas (ES), nine osteosarcomas (OS) and three dedifferentiated chondrosarcomas.
  
)

colAnn <- HeatmapAnnotation(df = ann_for_2_types,
                            which = "col",
                            col = colours,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1,"mm"))


#check order of Clinical and molecualr data
all.equal(colnames(t(mat_for_heatmap_ordered)),Ordered_pheno$Sample_ID)

#replace "." in name of the cell type with space 
colnames(mat_for_heatmap_ordered)<-gsub("."," ",fixed=TRUE,colnames(mat_for_heatmap_ordered))

ht4AB =Heatmap(t(mat_for_heatmap_ordered), name = "Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
             top_annotation =colAnn,show_column_dend = TRUE,show_row_dend = FALSE,show_column_names = FALSE)


#get expression values for genes that are part of PDl1 pathway and add them to heatmap 
mat_PDL1_related_genes_matrix<-normalized_counts[rownames(normalized_counts) %in% c("CTLA4","CD274","LAG3"),]

all.equal(colnames(mat_PDL1_related_genes_matrix),Ordered_pheno$Sample_ID)

mat_PDL1_related_genes_matrix_no_batch<-removeBatchEffect(mat_PDL1_related_genes_matrix, Ordered_pheno$batch,Ordered_pheno$`Abbreviation for Figures`)

ht4C =Heatmap(mat_PDL1_related_genes_matrix_no_batch,  name = "ICI genes",clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
               show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE)



##########. calculate TLS score

#From here: https://www.nature.com/articles/s41568-019-0144-6, 12 chemokine signature, and here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725838/
TLS_genes<-data.frame(c("CCL2","CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"))


colnames(TLS_genes)<-"TLS_genes"

##NOW calculate ssGSES TLS score for each sample using normalized counts 

TLS_genes<-data.frame(c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"))
colnames(TLS_genes)<-"TLS_genes"

TLS_genes_no_B_cell_related<-data.frame(c("CD1D", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"))
colnames(TLS_genes_no_B_cell_related)<-"TLS_genes"


##NOW calculate ssGSES TLS score for each sample using normalized counts 

#calculate TLS scores
geneSets_TLS<-as.list(TLS_genes)
geneSets_TLS$TLS_genes<-as.character(geneSets_TLS$`TLS signature`)
geneSets_TLS<-lapply(geneSets_TLS, function(z){ z[!is.na(z) & z != ""]})

gsva_TLS_gene<-gsva(as.matrix(normalized_counts),TLS_genes,method="gsva",kcdf="Gaussian" ,verbose=FALSE)

geneSets_TLS_no_B_cell_related<-as.list(TLS_genes_no_B_cell_related)
geneSets_TLS_no_B_cell_related$TLS_genes<-as.character(geneSets_TLS_no_B_cell_related$`TLS signature`)
geneSets_TLS_no_B_cell_related<-lapply(geneSets_TLS_no_B_cell_related, function(z){ z[!is.na(z) & z != ""]})

gsva_TLS_gene_no_B_cell_related<-gsva(as.matrix(normalized_counts),TLS_genes_no_B_cell_related,method="gsva",kcdf="Gaussian" ,verbose=FALSE)


#Function to plot regression
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

all.equal(colnames(normalized_counts),rownames(mat_for_heatmap_ordered),Ordered_pheno$CMO_ID,rownames(gsva_TLS_gene),rownames(gsva_TLS_gene_no_B_cell_related))

IKZF1_Bcell<-data.frame(normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),],mat_for_heatmap_ordered[,"B cell"],t(gsva_TLS_gene),t(gsva_TLS_gene_no_B_cell_related))
colnames(IKZF1_Bcell)<-c("IKZF1","B cell","TLS signature","TLS signature no B cell")

reg_IKZF1_B_cells<-ggplotRegression(lm(IKZF1_Bcell$`B cell`~ IKZF1_Bcell$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))
reg_IKZF1_B_cells_update<-reg_IKZF1_B_cells+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                             size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("B cell") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


reg_IKZF1_TLS<-ggplotRegression(lm(IKZF1_Bcell$IKZF1~ IKZF1_Bcell$`TLS signature`+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))

reg_IKZF1_TLS<-reg_IKZF1_TLS+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                              size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("TLS signature") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))



reg_IKZF1_TLS_noB<-ggplotRegression(lm(IKZF1_Bcell$IKZF1~ IKZF1_Bcell$`TLS signature no B cell`+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))

reg_IKZF1_TLS_noB<-reg_IKZF1_TLS_noB+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                      size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("TLS signature no B") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

##############################.  COMBINE ALL THE FIGURES


#COMBINE HEATMAPS vertically  
FIG4_ABC = ht4AB %v% ht4C
draw(FIG4_ABC)

dev.off()

############################################################################################################################

#             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR Immune hot/cold, Figure 4B

############################################################################################################################

#calculating KM For Immune hot/cold clusters


#readinf PFS daata
SARC028_PFS_data_trial<-fread("SARC028_PFS_data_trial.txt")


#Coding censored=0, 1=event occured

#importatnt nb of evens should be equal to nb of drops on the curve, each event is the drop on the curve. 
#and should be a cross , and each censpored should be tick so we have total of 3 ticks and they shuold be all in immune hot, 
# 

SARC028_PFS_data_trial$PFS_event<- ifelse(SARC028_PFS_data_trial$`PFS_Survival Event` %in% "Censored",0,1)


#join PFs with other clinical data
Ordered_pheno<-left_join(Ordered_pheno,SARC028_PFS_data_trial )



SARC028_baselinedata<-fread("SARC028_baselinedata.txt")

Ordered_pheno<-left_join(Ordered_pheno,SARC028_baselinedata )

Immune_hotSample_ID<-Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$Sample_ID
Immune_coldCSample_ID<-Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$Sample_ID


#test for age, BMI, or number of prior lines of therapy do not differ
t.test(Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_hotSample_ID,]$Age,Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_coldCSample_ID,]$Age)

t.test(Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_hotSample_ID,]$BMI,Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_coldCSample_ID,]$BMI)
t.test(Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_hotSample_ID,]$`Number of Prior Lines of Therapy`,Ordered_pheno[Ordered_pheno$Sample_ID %in% Immune_coldCSample_ID,]$`Number of Prior Lines of Therapy`)

Ordered_pheno$Responsebin<-ifelse(Ordered_pheno$Response %in% "CR/PR","responder","non-responder" )
#test for sex is also not signifincat 
fisher.test(table(Ordered_pheno$Sex,Ordered_pheno$Immune_subtypes))

#does not make much sense with only 3 responders
fisher.test(table(Ordered_pheno$Responsebin,Ordered_pheno$Immune_subtypes))

sfit <- survfit(Surv(PFS_months, PFS_event)~Immune_subtypes, data=Ordered_pheno)
ggsurvplot(sfit,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("#0000FF","red"),
                      legend.labs=c("Immune-cold","Immune-hot"),font.legend=14,
                      font.x=14,font.y=14,font.tickslab=14) 



model2 <- coxph( Surv(PFS_months, PFS_event) ~ Immune_subtypes  + `Abbreviation for Figures`, data = Ordered_pheno )
summary(model2)
ggforest(model2)




#Calaualte median PFS for immune hot/cold

median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$PFS_months)
median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$PFS_months)
############################################################################################################################



#for CD274
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotSample_ID][1,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCSample_ID][1,])

#for CTLA4
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotSample_ID][2,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCSample_ID][2,])

#for LAG3
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotSample_ID][3,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCSample_ID][3,])


##Calcauation of overall response rates (ORR) for immune hot and cold
#ORR=nb of responders
table(Ordered_pheno$Immune_subtypes,Ordered_pheno$Response)


#################################  correlation between IKZF1 and B-cells Supplementray Figure 10C ######################################
#Function to plot regression
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

all.equal(colnames(normalized_counts),rownames(mat_for_heatmap_ordered),Ordered_pheno$Sample_ID,rownames(gsva_TLS_gene))

IKZF1_Bcell<-data.frame(normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),],mat_for_heatmap_ordered[,"B cell"],t(gsva_TLS_gene))
colnames(IKZF1_Bcell)<-c("IKZF1","B cell","TLS signature")

reg_IKZF1_B_cells<-ggplotRegression(lm(IKZF1_Bcell$`B cell`~ IKZF1_Bcell$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))
reg_IKZF1_B_cells_update<-reg_IKZF1_B_cells+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                        size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("B cell") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


reg_IKZF1_TLS<-ggplotRegression(lm(IKZF1_Bcell$IKZF1~ IKZF1_Bcell$`TLS signature`+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))

reg_IKZF1_TLS<-reg_IKZF1_TLS+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                             size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("TLS signature") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

  

############################################################################################################################

#             GET Trasposable element counts from REDISCOVERTE, normalize and visualize TE counts

############################################################################################################################


##################################.   GET Transportable element counts from REDISCOVER TE pipeline and prepare for analysis  ##########################################
### get TEs expression only for LTR, LINE, SINE elements

#GEt only Intergenic TEs
INTERGENIC_TE_normalized_expression_1052_repeats<-read.table("INTERGENIC_TE_normalized_expression_1052_repeats")
colnames(INTERGENIC_TE_normalized_expression_1052_repeats)<-gsub(".","-",fixed=TRUE,colnames(INTERGENIC_TE_normalized_expression_1052_repeats))


#Match TEs with their familes 
TE_names_order<-data.frame(rownames(INTERGENIC_TE_normalized_expression_1052_repeats))
colnames(TE_names_order)<-"repName"
TE_family_annotation<-read.table("TE_family_annotation")
TE_names_order<-left_join(TE_names_order, TE_family_annotation)
TE_names_order$repClass<-gsub("LTR?", "LTR", fixed=TRUE, TE_names_order$repClass)
TE_names_order$repClass<-gsub("SINE?", "SINE", fixed=TRUE, TE_names_order$repClass)

#seelct for samples that we have information for 
INTERGENIC_TE_normalized_expression_1052_repeats<-INTERGENIC_TE_normalized_expression_1052_repeats[,colnames(INTERGENIC_TE_normalized_expression_1052_repeats) %in% Ordered_pheno$Sample_ID]

#### Filter gene TE matrix filterByExpr
#TEs are rows, samples are columns, use expression per group as a minimum 

#select those that are baseline 
INTERGENIC_TE_normalized_expression_1052_repeats<-INTERGENIC_TE_normalized_expression_1052_repeats[,colnames(INTERGENIC_TE_normalized_expression_1052_repeats) %in% Ordered_pheno$Sample_ID]

#check the order of samples
all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats), Ordered_pheno$Sample_ID)



keep.exprs_TE <- filterByExpr(INTERGENIC_TE_normalized_expression_1052_repeats, group=Ordered_pheno$`Abbreviation for Figures`)
INTERGENIC_TE_normalized_expression_1052_repeats_filtered <- INTERGENIC_TE_normalized_expression_1052_repeats [which(keep.exprs_TE==TRUE),]
dim(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)

all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_filtered), Ordered_pheno$Sample_ID)

# read distribution across TE families;
rowmeans_no_filters<-rowMeans(INTERGENIC_TE_normalized_expression_1052_repeats)
Hist_TE<-hist(rowmeans_no_filters)

median(rowmeans_no_filters)


rowmeans_filters<-rowMeans(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)
hist(rowmeans_filters,xlim =c(0,8000))
median(rowmeans_filters)

dim(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)


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
all.equal(Ordered_pheno$Sample_ID,colnames(INTERGENIC_TE_normalized_expression_1052_repeats))

##################### PLOT INTERGENIC NORMALIZED EXPRESSION OVER TWO CLSUTERS

INTERGENIC_TE_normalized_expression_1052_repeats_no_batch<-removeBatchEffect(INTERGENIC_TE_normalized_expression_1052_repeats_Z_sore_transformed, Ordered_pheno$`Abbreviation for Figures`)

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


###Create annotations and plot

#Annotation for columns
ann_for_2_types <- data.frame(Ordered_pheno[,c("Response","Abbreviation for Figures")]) 
rownames(ann_for_2_types)<-Ordered_pheno$"Sample_ID"
colnames(ann_for_2_types) <- c("Response","Subtype") 



colours <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  'Subtype'=c ( "DDLPS"= "#DEABCA","UPS"= "#82D992", "LS"= "#F4A460","OS"= "#8881D2" , "CS"="green","ES"="purple","SS"="magenta")
  
  
)


colAnn2 <- HeatmapAnnotation(df = ann_for_2_types,
                             which = "col",
                             col = colours,
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

Heatmap(TE_intergenic_for_heatmap, name = "mcp_counter_Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
        top_annotation =colAnn2,show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE,cluster_columns =TRUE,show_row_names = FALSE,right_annotation = rowAnn,cluster_rows = FALSE,row_order=as.numeric(rownames(TE_intergenic_for_heatmap_row_anno[order(TE_intergenic_for_heatmap_row_anno$repClass),])))



sessionInfo()

# to create a report from this script just open Rstudio and render the script using rmarkdown, note: the script and data must be all in the same folder
#The commnad to render is rmarkdown::render("Script_1.R")

