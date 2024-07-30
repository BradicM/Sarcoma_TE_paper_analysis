#Author Martina Bradic
#This script performs a comprehensive analysis of our clinical raw RNA sequencing counts obtained from the REDISCOVERTE
#pipeline run. The script first normalizes and filters the counts, and then performs immune deconvolution using the MCP 
#counter. The deconvoluted immune cell proportions are then clustered using FactoMineR to establish two major
#clusters: immune hot and cold.
#The script generates several figures, including a heatmap of the immune clusters and relevant clinical correlates 
#(Figure 1A,B,C,D), a Kaplan-Meier analysis and plot for Immune hot/cold groups and PFS (months) (Figure 1E),
#factor miner clustering of deconvoluted data (Supplemental Figure 1), a Cox Hazardous model analysis and plot
#(Supplemental Figure 2), and a matching of our immune hot and cold clusters with previously published 
#Petitprez et.al clustering of immune cells (Supplemental Figure 3). Additionally, the script generates a 
#heatmap of normalized intergenic TE counts with clinical correlates (Supplemental Figure 5).
#The script also includes several supplementary figures, including a correlation between IKZF1
#and B-cells (Supplemental Figure 6A).

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

############ GET gene expression
gene_expression_mtrx<-read.table("gene_expression_mtrx")
colnames(gene_expression_mtrx)<-gsub(".","-",fixed=TRUE,colnames(gene_expression_mtrx))

#get clinical data
Clinical_phenotypes<-fread("Clinical_phenotypes.txt")


#order clinical file based on gene expresion 
Ordered_pheno<-data.frame(colnames(gene_expression_mtrx))
colnames(Ordered_pheno)<-"CMO_ID"
Ordered_pheno<-left_join(Ordered_pheno,Clinical_phenotypes)

#perform filtering
keep.exprs <- filterByExpr(gene_expression_mtrx, group=Ordered_pheno$`Abbreviation for Figures`)
gene_expression_mtrx_filtered <- gene_expression_mtrx[which(keep.exprs==TRUE),]
dim(gene_expression_mtrx_filtered)

#    after filtering 
#30699    67

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




########################  PLOTTING HEATMAPS ##########################################
#we remove batch before
mat_for_heatmap<-data.frame(t(mat_immuno2_no_batch))

mat_for_heatmap$CMO_ID<-rownames(mat_for_heatmap)
mat_for_heatmap_ordered<-data.frame(Ordered_pheno$CMO_ID)
colnames(mat_for_heatmap_ordered)<-"CMO_ID"
mat_for_heatmap_ordered<-left_join(mat_for_heatmap_ordered, mat_for_heatmap)
rownames(mat_for_heatmap_ordered)<-mat_for_heatmap_ordered$CMO_ID
mat_for_heatmap_ordered<-mat_for_heatmap_ordered[,c(-1)]

ann_for_2_types <- data.frame(Ordered_pheno[,c("Response","Abbreviation for Figures")]) 
rownames(ann_for_2_types)<-Ordered_pheno$CMO_ID
colnames(ann_for_2_types) <- c("Response","Subtype") 

colours <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  'Subtype'=c ("ANGS"="#DD67C0" ,"CHS"="#A848E0","DDLS"= "#DEABCA","LMS"= "#F4A460","OS"= "#8881D2" ,"SARCNOS"="#8DB6D1","UPS"= "#82D992","EHE"="#FFFF00","LPS"="#0000FF","MFS"="#B22222","Other"="#FF0000","SBRC"="#85E0D6","ASPS"="grey")
  
  
)

colAnn <- HeatmapAnnotation(df = ann_for_2_types,
                            which = "col",
                            col = colours,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1,"mm"))


#check order of Clinical and molecualr data
all.equal(colnames(t(mat_for_heatmap_ordered)),Ordered_pheno$CMO_ID)

#replace "." in name of the cell type with space 
colnames(mat_for_heatmap_ordered)<-gsub("."," ",fixed=TRUE,colnames(mat_for_heatmap_ordered))

ht1AB =Heatmap(t(mat_for_heatmap_ordered), name = "Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
             top_annotation =colAnn,show_column_dend = TRUE,show_row_dend = FALSE,show_column_names = FALSE)


#get expression values for genes that are part of PDl1 pathway and add them to heatmap 
mat_PDL1_related_genes_matrix<-normalized_counts[rownames(normalized_counts) %in% c("CTLA4","CD274","LAG3"),]

all.equal(colnames(mat_PDL1_related_genes_matrix),Ordered_pheno$CMO_ID)

mat_PDL1_related_genes_matrix_no_batch<-removeBatchEffect(mat_PDL1_related_genes_matrix, Ordered_pheno$batch,Ordered_pheno$`Abbreviation for Figures`)

ht1C =Heatmap(mat_PDL1_related_genes_matrix_no_batch,  name = "ICI genes",clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_2_types),column_split =Ordered_pheno$Immune_subtypes,
               show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE)


##############################.  COMBINE ALL THE FIGURES


FIG1_ABC = ht1AB %v% ht1C 
draw(FIG1_ABC)


############################################################################################################################

#             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR Immune hot/cold, Figure 1E and Supplemental Fig 2

############################################################################################################################

#calculating KM For Immune hot/cold clusters

sfit <- survfit(Surv(PFS_months, `Alive/dead`)~Immune_subtypes, data=Ordered_pheno)
ggsurvplot(sfit,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("#0000FF","red"),
                      legend.labs=c("Immune-cold","Immune-hot"),font.legend=14,
                      font.x=14,font.y=14,font.tickslab=14) 



#Supplemental Fig 2 coxph anlysis for Immune subtype
#relevel to make UPS as a reference
Ordered_pheno_relevel<-Ordered_pheno
Ordered_pheno_relevel$`Abbreviation for Figures`<- as.factor(Ordered_pheno_relevel$`Abbreviation for Figures` )
Ordered_pheno_relevel$`Abbreviation for Figures`<-relevel(Ordered_pheno_relevel$`Abbreviation for Figures`, ref = "OS")

model2 <- coxph( Surv(PFS_months, `Alive/dead`) ~ Immune_subtypes  + `Abbreviation for Figures`, data = Ordered_pheno_relevel )
ggforest(model2)



#Calaualte median PFS for immune hot/cold

median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$PFS_months)
median(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$PFS_months)
############################################################################################################################

#             FISHER EXACT ANALYSIS COMPARING CLINICAL PROTOCOLS FOR RESPONDERS AND IMMUNE HOT/COLD
# These analysis is summarized in Supplemental table 3 and in the text
############################################################################################################################

#Perform fisher exact test to determine if there is differences in number of responders/non-reseponders in hot/cold immune types

fisher.test(table(Ordered_pheno$Immune_subtypes,Ordered_pheno$Response2))

fisher.test(table(Ordered_pheno[Ordered_pheno$Response %in% c("CR/PR","PD"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Response %in% c("CR/PR","PD"),]$Response))
fisher.test(table(Ordered_pheno[Ordered_pheno$Response %in% c("CR/PR","SD"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Response %in% c("CR/PR","SD"),]$Response))
fisher.test(table(Ordered_pheno[Ordered_pheno$Response %in% c("PD","SD"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Response %in% c("PD","SD"),]$Response))



fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-366"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-366"),]$Protocol))

fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("17-508","17-366"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Protocol %in% c("17-508","17-366"),]$Protocol))

fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-508"),]$Immune_subtypes,Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-508"),]$Protocol))


fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-366"),]$Response2 ,Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-366"),]$Protocol))

fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("17-508","17-366"),]$Response2,Ordered_pheno[Ordered_pheno$Protocol %in% c("17-508","17-366"),]$Protocol))

fisher.test(table(Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-508"),]$Response2,Ordered_pheno[Ordered_pheno$Protocol %in% c("16-1534","17-508"),]$Protocol))


table(Ordered_pheno$Immune_subtypes,Ordered_pheno$Protocol)


##Also  t.test for expression levels of immune checkpoint-related genes  from Figure 1C

Immune_hotCMOIDs<-Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$CMO_ID
Immune_coldCMOIDs<-Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$CMO_ID

#for CD274
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotCMOIDs][1,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCMOIDs][1,])

#for CTLA4
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotCMOIDs][2,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCMOIDs][2,])

#for LAG3
t.test(mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_hotCMOIDs][3,],mat_PDL1_related_genes_matrix_no_batch[,colnames(mat_PDL1_related_genes_matrix_no_batch) %in% Immune_coldCMOIDs][3,])


##Calcauation of overall response rates (ORR) for immune hot and cold
#ORR=nb of responders
table(Ordered_pheno$Immune_subtypes,Ordered_pheno$Response2)


#################################  Supplemental Figure 6A, correlation between IKZF1 and TEs with TLS, and B cells ######################################

##########. calculate TLS score

TLS_genes<-data.frame(c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"))
colnames(TLS_genes)<-"TLS_genes"

TLS_genes_no_B_cell_related<-data.frame(c("CD1D", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"))
colnames(TLS_genes_no_B_cell_related)<-"TLS_genes"


##NOW calculate ssGSES TLS score for each sample using normalized counts 

#calculate TLS scores

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

reg_IKZF1_B_cells<-ggplotRegression(lm(IKZF1_Bcell$`B cell`~ IKZF1_Bcell$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))
reg_IKZF1_B_cells_update<-reg_IKZF1_B_cells+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                        size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("B cell") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


reg_IKZF1_TLS<-ggplotRegression(lm(IKZF1_Bcell$IKZF1~ IKZF1_Bcell$`TLS signature`+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))

reg_IKZF1_TLS<-reg_IKZF1_TLS+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                             size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("TLS signature") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

  

reg_IKZF1_TLS_noB<-ggplotRegression(lm(IKZF1_Bcell$IKZF1~ IKZF1_Bcell$`TLS signature no B cell`+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))

reg_IKZF1_TLS_noB<-reg_IKZF1_TLS_noB+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                              size=14),axis.text.y = element_text(face="bold", size=14)) +
  xlab('IKZF1') + ylab("TLS signature no B") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

############################################################################################################################

#         Comparison of this paper's Immune types with published work from Petitpreze et.al.,
## This portion of code performs analysis and plots heatmap of matchcing clusters Supplemental Figure 3

############################################################################################################################

#read in centroid table that they have send from Petitprez paper, and match Immune hot and cold types
# using euclidian distance to determine which samples belong to which published (SIC)

Petitprez_cluster_centroids<-read.table("NanostringCentroids_from_Petiprez_publication.txt", sep="\t")
rownames(Petitprez_cluster_centroids)<-Petitprez_cluster_centroids$V1
Petitprez_cluster_centroids<-Petitprez_cluster_centroids[-c(1),-c(1)]
colnames(Petitprez_cluster_centroids)<-c("A","B","C","D","E")

For_euclidian_distance<-mat_for_heatmap_ordered[,c("T cell","cytotoxicity score","B cell" ,"Endothelial cell")]

For_euclidian_distance<-t(For_euclidian_distance)


euclidean <- function(a, b) sqrt(sum((a - b)^2))


# now calculate euclidean distance between
SIC_A<-list()
for (i in 1:ncol(For_euclidian_distance)){
  SIC_A[i]<-euclidean(For_euclidian_distance[,i],as.numeric(as.character(Petitprez_cluster_centroids[,1])))       
}

SIC_A<-data.matrix(SIC_A)

SIC_B<-list()
for (i in 1:ncol(For_euclidian_distance)){
  SIC_B[i]<-euclidean(For_euclidian_distance[,i],as.numeric(as.character(Petitprez_cluster_centroids[,2])))       
}

SIC_B<-data.matrix(SIC_B)

SIC_C<-list()
for (i in 1:ncol(For_euclidian_distance)){
  SIC_C[i]<-euclidean(For_euclidian_distance[,i],as.numeric(as.character(Petitprez_cluster_centroids[,3])))      
}

SIC_C<-data.matrix(SIC_C)


SIC_D<-list()
for (i in 1:ncol(For_euclidian_distance)){
  SIC_D[i]<-euclidean(For_euclidian_distance[,i],as.numeric(as.character(Petitprez_cluster_centroids[,4])))      
}

SIC_D<-data.matrix(SIC_D)

SIC_E<-list()
for (i in 1:ncol(For_euclidian_distance)){
  SIC_E[i]<-euclidean(For_euclidian_distance[,i],as.numeric(as.character(Petitprez_cluster_centroids[,5])))      
}

SIC_E<-data.matrix(SIC_E)


SIC_assignemnt_results<-cbind(SIC_A,SIC_B,SIC_C,SIC_D, SIC_E)

colnames(SIC_assignemnt_results)<-c("A","B","C","D","E")
rownames(SIC_assignemnt_results)<-colnames(For_euclidian_distance)
SIC_assignemnt_results<-data.frame(SIC_assignemnt_results)
SIC_assignemnt_results$min_col = colnames(SIC_assignemnt_results)[apply(SIC_assignemnt_results, 1, which.min)]

# Assign SIC assigment to phenotye file for heatmap
SIC_assignemnt_results$CMO_ID<-rownames(SIC_assignemnt_results)
Ordered_pheno_SIC<-left_join(Ordered_pheno,SIC_assignemnt_results[,c("CMO_ID","min_col")], by=c("CMO_ID"))
colnames(Ordered_pheno_SIC)[20]<-"SIC clusters"



ann3 <- data.frame(Ordered_pheno_SIC[,c("Response","Abbreviation for Figures", "Immune_subtypes","SIC clusters")]) 

rownames(ann3)<-Ordered_pheno_SIC$CMO_ID
colnames(ann3) <- c("Response","Subtype","Immune types","SICs") 


colours3 <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  'Subtype'=c ("ANGS"="#DD67C0" ,"CHS"="#A848E0","DDLS"= "#DEABCA","LMS"= "#F4A460","OS"= "#8881D2" ,"SARCNOS"="#8DB6D1","UPS"= "#82D992","EHE"="#FFFF00","LPS"="#0000FF","MFS"="#B22222","Other"="#FF0000","SBRC"="#85E0D6","ASPS"="grey"),
  'Immune types'=c("Immune hot"="red","Immune cold"="blue")
  
  
)

col3 <- HeatmapAnnotation(df = ann3,
                          which = "col",
                          col = colours3,
                          annotation_width = unit(c(1, 4), "cm"),
                          gap = unit(1,"mm"))



Heatmap(t(mat_for_heatmap_ordered), name = "mcp_counter_Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann3),column_split =Ordered_pheno_SIC$`SIC clusters`,
        top_annotation =col3,show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE,row_order=c(10,11,3,9,1,6,5,7,8,4,2), cluster_columns = FALSE)


#count how many immune hot/cold are in the published clusters from Petitprez paper
table(Ordered_pheno_SIC$Immune_subtypes, Ordered_pheno_SIC$`SIC clusters`)



############################################################################################################################

#             GET Trasposable element counts from REDISCOVERTE, normalize and visualize as heatmap (Supplemental Fig 5)

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


#### Filter gene TE matrix filterByExpr
#TEs are rows, samples are columns, use expression per group as a minimum 

#check the order of samples
all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats), Ordered_pheno$CMO_ID)

keep.exprs_TE <- filterByExpr(INTERGENIC_TE_normalized_expression_1052_repeats, group=Ordered_pheno$`Abbreviation for Figures`)
INTERGENIC_TE_normalized_expression_1052_repeats_filtered <- INTERGENIC_TE_normalized_expression_1052_repeats[which(keep.exprs_TE==TRUE),]
dim(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)

all.equal(colnames(INTERGENIC_TE_normalized_expression_1052_repeats_filtered), Ordered_pheno$CMO_ID)

# read distribution across TE families; Supplemental Figure 7
rowmeans_no_filters<-rowMeans(INTERGENIC_TE_normalized_expression_1052_repeats)
hist(rowmeans_no_filters,xaxt='n')
axis(side=1,at=seq(0,8000,500),labels=seq(0,8000,500))

median(rowmeans_no_filters)

rowmeans_filters<-rowMeans(INTERGENIC_TE_normalized_expression_1052_repeats_filtered)

hist(rowmeans_filters, xaxt='n')
#axis(side=1, at=seq(0,500, 1000,1500,2000), labels=seq(0,500, 1000,1500,2000))

axis(side=1,at=seq(0,8000,500),labels=seq(0,8000,500))
median(rowmeans_filters)



###There are 1002 TEs that are left here to be analysed

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
ann_for_2_types <- data.frame(Ordered_pheno[,c("Response","Abbreviation for Figures")]) 
rownames(ann_for_2_types)<-Ordered_pheno$"CMO_ID"
colnames(ann_for_2_types) <- c("Response","Subtype") 



colours <- list(
  'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
  'Subtype'=c ("ANGS"="#DD67C0" ,"CHS"="#A848E0","DDLS"= "#DEABCA","LMS"= "#F4A460","OS"= "#8881D2" ,"SARCNOS"="#8DB6D1","UPS"= "#82D992","EHE"="#FFFF00","LPS"="#0000FF","MFS"="#B22222","Other"="#FF0000","SBRC"="#85E0D6","ASPS"="grey")
  
  
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
#The command to render is rmarkdown::render("Script_1.R")

