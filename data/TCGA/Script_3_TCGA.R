#Author Martina Bradic
#This script requires that you first run Figure 1 script
#This script gets gene lists for immune pathways from published data, calculates TE-score for significant TE familes, 
#and performs and partical correlation and plotting of Figure 9A, 
# And it also performs Kaplan Meier plotting and analysis for Figure 9B (for IKZF1 and TE score)
library(data.table)
library(survminer)
library(survival)
library(dplyr)
library(dotwhisker)
library(GSVA)
library(ppcor)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(corrplot)
library(tidyr)

############################################################################################################################

#             GET immune pathway and C-gas genes for GSVA calculation, and calcualte TE score

############################################################################################################################
 


### calcaute ssGSVA for candidate gene lists from Nature communication paper https://www.nature.com/articles/s41467-019-13035-2
##Get gene list

Gene_lists_for_ddGSEA<-fread("Gene_lists_for_ddGSEA.txt",header=TRUE)


Gene_lists_for_ddGSEA<-Gene_lists_for_ddGSEA %>% 
  mutate(Genes=strsplit(Genes, ",")) %>% 
  unnest(Genes)

Gene_list<-split.data.frame(Gene_lists_for_ddGSEA[-1],Gene_lists_for_ddGSEA$`Gene Signature`)
geneSets<-lapply(Gene_list,"[[","Genes")


##NOW calculate ssGSES for each sample using normalized counts 
gsva_24_gene_lists<-gsva(as.matrix(normalized_counts),geneSets,method="gsva",kcdf="Gaussian" ,verbose=FALSE)

#### Associations with that list

all.equal(colnames(gsva_24_gene_lists),Ordered_pheno$patient)

Gene_list_24_immune_clustering<-vector(mode = "list", length = ncol(t(gsva_24_gene_lists)))
for (i in 1:ncol(t(gsva_24_gene_lists))) {
  fit_Test<-lm( t(gsva_24_gene_lists)[,i] ~ as.numeric(as.factor(Ordered_pheno$Immune_subtypes)) + Ordered_pheno$batch +  Ordered_pheno$`Abbreviation for Figures`)  # 
  Gene_list_24_immune_clustering[[i]]<-anova(fit_Test)$'Pr(>F)'[1]
}

Gene_list_24_immune_clustering_df<-data.frame(do.call("rbind",Gene_list_24_immune_clustering) )
rownames(Gene_list_24_immune_clustering_df)<-rownames(gsva_24_gene_lists)
colnames(Gene_list_24_immune_clustering_df)<-"Pvalue"
Gene_list_24_immune_clustering_df$adjPvalue<-p.adjust(Gene_list_24_immune_clustering_df$Pvalue, method ="BH",n =length(Gene_list_24_immune_clustering_df$Pvalue))

#get only significant once
Gene_list_24_immune_clustering_df[Gene_list_24_immune_clustering_df$adjPvalue<0.05,]

gsva_24_gene_lists2<-data.frame(gsva_24_gene_lists)

#convert to z-score
mat_gsva_24_gene_list<-scale(gsva_24_gene_lists2)

colnames(mat_gsva_24_gene_list)<-gsub(".","-",fixed=TRUE,colnames(mat_gsva_24_gene_list))
#remove batch (batch + histology)
#are orders of the samples the same, yes
all.equal(colnames(mat_gsva_24_gene_list),Ordered_pheno$patient)

mat_gsva_24_gene_list<-removeBatchEffect(mat_gsva_24_gene_list, Ordered_pheno$batch)

### Get c gas genes


C_GAS_gene_list_KEGG<-fread("C_GAS_gene_list_KEGG")
C_GAS_gene_list_KEGG<-C_GAS_gene_list_KEGG[-c(1),]
C_GAS_gene_list_KEGG<-data.frame(C_GAS_gene_list_KEGG)
colnames(C_GAS_gene_list_KEGG)<-"genes"


##Calculate ssGSES for each sample using normalized counts 
gsva_result_C_gas<-gsva(as.matrix(normalized_counts),list(C_GAS_gene_list_KEGG$genes),method="gsva",kcdf="Gaussian",verbose=FALSE )

gsva_result_C_gas_mat<-data.frame(t(gsva_result_C_gas))
gsva_result_C_gas_mat_Z_sore_transformed<-t(gsva_result_C_gas_mat %>% mutate_at(colnames(gsva_result_C_gas_mat), scale))
colnames(gsva_result_C_gas_mat_Z_sore_transformed)<-colnames(gsva_result_C_gas)


##Calculate correlations between TE features (Supplemental Figure 7A) and TE score

TE_features<-c(rownames(TE_and_Cellularity_model_features$glmnet_results_df)[-c(1)])
TE_features<-data.frame(TE_features)



TE_genes_df<-as.data.frame(TE_matrix[,colnames(TE_matrix) %in% as.character(TE_features$TE_features)])
M_TE_corr <-cor(TE_genes_df)


#correlation between significant TEs
corrplot(M_TE_corr, type="upper",method = 'circle', order = 'AOE',tl.col = "black",tl.cex = 0.8,
         col=colorRampPalette(c("blue","white","red"))(200))


#calculate TE scores
geneSets<-as.list(TE_features)
geneSets$TE_features<-as.character(geneSets$TE_features)
geneSets<-lapply(geneSets, function(z){ z[!is.na(z) & z != ""]})

TE_genes_matrix<-as.matrix(TE_matrix[,colnames(TE_matrix) %in% as.character(TE_features$TE_features)])

gsva_for_TE__genes <-gsva(t(TE_genes_matrix),geneSets,method="zscore", kcdf="Gaussian",verbose=FALSE)

#check for same order
all.equal(colnames(gsva_for_TE__genes),colnames(mat_gsva_24_gene_list))

#get all the data into the same matrix for downstream calculation 
mat3<-data.frame(rbind(mat_gsva_24_gene_list,gsva_result_C_gas_mat$t.gsva_result_C_gas.,normalized_counts_all_EPI_genes["IKZF1",],gsva_for_TE__genes,normalized_counts[c("CTLA4","CD274"),]))
rownames(mat3)[c(25,26,27)]<-c("C-gas_pathway","IKZF1","TE score")
colnames(mat3)<-gsub(".","-",fixed=TRUE,colnames(mat3))
mat3_transpose<-data.frame(t(mat3))

all.equal(rownames(mat3_transpose),Ordered_pheno$patient)



############################################################################################################################

#           Calculate partial correlation between pathways and CD274,TE score and IKZF1 and plot Supplemental figure 9A 

############################################################################################################################


#Pcor test for all the pairwise comparison 
estimate<-sapply(1:(ncol(mat3_transpose)), function(x) sapply(1:(ncol(mat3_transpose)), function(y) {
  if (x == y) 1
  else pcor.test(mat3_transpose[,x], mat3_transpose[,y],method = c("pearson"), list(as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`)),as.numeric(as.factor(Ordered_pheno$batch))))$estimate
}))

colnames(estimate)<-colnames(mat3_transpose)
rownames(estimate)<-colnames(mat3_transpose)

pvalues<-sapply(1:(ncol(mat3_transpose)), function(x) sapply(1:(ncol(mat3_transpose)), function(y) {
  if (x == y) 1
  else pcor.test(mat3_transpose[,x], mat3_transpose[,y],method = c("pearson"), list(as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`)),as.numeric(as.factor(Ordered_pheno$batch))))$p.value
}))

colnames(pvalues)<-colnames(mat3_transpose)
rownames(pvalues)<-colnames(mat3_transpose)

#now  select only unique columns/rows
estimate<-estimate[-c(26:29),c(26:29)]
pvalues<-pvalues[-c(26:29),c(26:29)]

#reformat data
melted_cormat <- melt(t(estimate))

melted_p <- melt(t(pvalues))
melted_p$BH_adjusted<-p.adjust(melted_p$value,method="BH")
melted_p$stars <- cut(melted_p$BH_adjusted, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

#join p value and correlation 
melted_cormat<-data.frame(melted_cormat,melted_p$stars )

#get pathway annotation table
Pathway_annotations_for_24_pathways<-fread("Pathway_annotations_for_24_pathways.txt")
melted_cormat<-left_join(melted_cormat,Pathway_annotations_for_24_pathways )

#remove CTLA4
melted_cormat<-melted_cormat[melted_cormat$Var1 %nin% "CTLA4",]


# FIGURE 10A plot ; Heatmap of partial correlations; note, the figure in paper in transpose orientation 

colnames(melted_cormat)<-c("Var1", "Pathway_name" ,"value", "melted_p.stars")

print(melted_cormat)



melted_cormat$ord.x <- factor(melted_cormat$Pathway_name, ordered=TRUE, levels = c("WNT.target","Pan.F.TBRS","P53.signalling.pathway",
                                                                                   "Nucleotide.excision.repair","NHEJ","Mismatch.repair","Homologous.recombination", "FGFR3.related.genees","Fanconi.anemia",
                                                                                   "EMT", "DNA.replication.dependent.histones","DNA.replication","DNA.damage.repair","Cell.cycle.regulators" ,"Cell.cycle" ,
                                                                                   "Angiogenesis","Type.II.IFN.Response","Type.I.IFN.Response","TNFalpha.Response","NFkB.response","Immune.checkpoint",
                                                                                   "IL1beta.Response","CD8.T.effector","C.gas_pathway","Antigen.processing.machinery"))


SupplementalFig9A<-ggplot(data = melted_cormat, aes(ord.x, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson partial correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+ 
  geom_text(aes(label=melted_p.stars), color="black", size=5)

SupplementalFig9A + coord_flip()
  

############################################################################################################################

#           Calculate correlation and  Plot Figure 9B 

############################################################################################################################

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



#calculate and plot regression for CD274 vs Te score, and CD274 vs IKZF1
#lm for CD274 and TEscore
summary(lm(mat3_transpose$CD274~ mat3_transpose$TE.score+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))

#plot the correlation 
reg1_CD274_TE<-ggplotRegression(lm(mat3_transpose$CD274~ mat3_transpose$TE.score+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))

reg1_CD274_TE_update<-reg1_CD274_TE+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                     size=14),axis.text.y = element_text(face="bold",
                                                                                                                                                                                         size=14)) +
  xlab("TE score") + ylab("Normalized expression CD274") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

#lm for CD274 and IKZF1
summary(lm(mat3_transpose$CD274~ mat3_transpose$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))
reg2_CD274_IKZF1<-ggplotRegression(lm(mat3_transpose$CD274~ mat3_transpose$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))
reg2_CD274_IKZF1_update<-reg2_CD274_IKZF1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                            size=14),axis.text.y = element_text(face="bold",                                                                                                                                                                      size=14)) +
  xlab("Normalized expression IKZF1") + ylab("Normalized expression CD274") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


ggarrange(reg1_CD274_TE_update, reg2_CD274_IKZF1_update,
          ncol = 1, nrow = 2, align = "v")


############################################################################################################################

#             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR IKZF1; Figure 4E

############################################################################################################################

#################. Calculating KM for TE score; ; Figure 4E
#check order

all.equal(rownames(t(gsva_for_TE__genes)),Ordered_pheno$patient)

For_KM_TEscore<-data.frame(cbind(t(gsva_for_TE__genes),Ordered_pheno))

TE_features.cut <- surv_cutpoint(
  data=For_KM_TEscore,
  time = "OS_months",
  event = "OS_status",
  variables = "TE_features"
)



For_KM_TEscore$TE_features_levels<-ifelse(For_KM_TEscore$TE_features> summary(TE_features.cut)[1]$cutpoint, "high TE score","low TE score")

sfit_TE <- survfit(Surv(OS_months, OS_status)~TE_features_levels , data=For_KM_TEscore)


ggsurvplot(sfit_TE,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
           legend.labs=c("high TE score","low TE score"),font.legend=14,
           font.x=14,font.y=14,font.tickslab=14) 


fit_TEscore <- coxph(Surv(OS_months, OS_status)~TE_features_levels  + Abbreviation.for.Figures + pathologic.tumor.size, data=For_KM_TEscore)
summary(fit_TEscore)



######################## Correlation between IKZF1 expression and B-cell proportions; this is not featured  in n the manuscript but it 
## also replicates the observation that IKZF1 and B-cells are significantly positively correlated
all.equal(rownames(mat_for_heatmap_ordered),Ordered_pheno$patient,colnames(normalized_counts))

mat_for_heatmap_ordered<-data.frame(mat_for_heatmap_ordered)
reg1<-ggplotRegression(lm(mat_for_heatmap_ordered$B.cell ~ normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),]+ as.numeric(as.factor(Ordered_pheno$batch))))

reg1_update<-reg1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                   size=14),axis.text.y = element_text(face="bold",
                                                                                                                                                                       size=14)) +
  xlab("IKZF1") + ylab("B cell") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

