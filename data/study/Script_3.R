#Author Martina Bradic
#This script requires that you first run Script_1.R and Script_2.R.
#This script retrieves gene lists for immune pathways from published data and calculates TE-scores for 
#significant TE families. It then performs  partial correlation and plotting of Figure 3A, 
#correlation and plotting of Figure 3B,  and Kaplan-Meier plotting and analysis for Figure 3C (for IKZF1 and TE score).

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
library(bnlearn)
library(dotwhisker)
library(jtools)
############################################################################################################################

#             GET immune pathway and C-gas genes for GSVA calculation

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

all.equal(colnames(gsva_24_gene_lists),Ordered_pheno$CMO_ID,colnames(gsva_TLS_gene))

Gene_list_24_immune_clustering<-vector(mode = "list", length = ncol(t(gsva_24_gene_lists)))
for (i in 1:ncol(t(gsva_24_gene_lists))) {
  fit_Test<-lm( t(gsva_24_gene_lists)[,i] ~ as.numeric(as.factor(Ordered_pheno$Immune_subtypes)) + Ordered_pheno$batch +  Ordered_pheno$`Abbreviation for Figures`)  # 
  Gene_list_24_immune_clustering[[i]]<-anova(fit_Test)$'Pr(>F)'[1]
}

Gene_list_24_immune_clustering_df<-data.frame(do.call("rbind",Gene_list_24_immune_clustering) )
rownames(Gene_list_24_immune_clustering_df)<-rownames(gsva_24_gene_lists)
colnames(Gene_list_24_immune_clustering_df)<-"Pvalue"
Gene_list_24_immune_clustering_df$adjPvalue<-p.adjust(Gene_list_24_immune_clustering_df$Pvalue, method ="BH",n =length(Gene_list_24_immune_clustering_df$Pvalue))

#10 out of 24 pathways are significant  by lm with batch and histology as a covariats
Gene_list_24_immune_clustering_df[Gene_list_24_immune_clustering_df$adjPvalue<0.05,]

gsva_24_gene_lists2<-data.frame(gsva_24_gene_lists)

#convert to z-score
mat_gsva_24_gene_list<-scale(gsva_24_gene_lists2)

colnames(mat_gsva_24_gene_list)<-gsub(".","-",fixed=TRUE,colnames(mat_gsva_24_gene_list))
#remove batch (batch + histology)
#are orders of the samples the same, yes
all.equal(colnames(mat_gsva_24_gene_list),Ordered_pheno$CMO_ID)

mat_gsva_24_gene_list<-removeBatchEffect(mat_gsva_24_gene_list, Ordered_pheno$batch,Ordered_pheno$`Abbreviation for Figures`)

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




############################################################################################################################

#           Calculate partial correlation between pathways and CD274,TE score and IKZF1 and plot Figure 3A 

############################################################################################################################
#prepare the data for the calcualation 

all.equal(colnames(gsva_for_TE__genes),colnames(mat_gsva_24_gene_list))

#get all the data into the same matrix for downstream calculation 
mat3<-data.frame(rbind(mat_gsva_24_gene_list,gsva_result_C_gas_mat$t.gsva_result_C_gas.,normalized_counts_all_EPI_genes["IKZF1",],gsva_for_TE__genes,normalized_counts[c("CTLA4","CD274"),]))
rownames(mat3)[c(25,26,27)]<-c("C-gas_pathway","IKZF1","TE score")
colnames(mat3)<-gsub(".","-",fixed=TRUE,colnames(mat3))
mat3_transpose<-data.frame(t(mat3))

all.equal(rownames(mat3_transpose),Ordered_pheno$CMO_ID)

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

#show partial correlation values
print(estimate)

#show partial p values
print(pvalues)
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


# FIGURE 3A plot ; Heatmap of partial correlations; note, the figure in paper in transpose orientation 

colnames(melted_cormat)<-c("Var1", "Pathway_name" ,"value", "melted_p.stars")

print(melted_cormat)



melted_cormat$ord.x <- factor(melted_cormat$Pathway_name, ordered=TRUE, levels = c("WNT.target","Pan.F.TBRS","P53.signalling.pathway",
                                                                                   "Nucleotide.excision.repair","NHEJ","Mismatch.repair","Homologous.recombination", "FGFR3.related.genees","Fanconi.anemia",
                                                                                   "EMT", "DNA.replication.dependent.histones","DNA.replication","DNA.damage.repair","Cell.cycle.regulators" ,"Cell.cycle" ,
                                                                                   "Angiogenesis","Type.II.IFN.Response","Type.I.IFN.Response","TNFalpha.Response","NFkB.response","Immune.checkpoint",
                                                                                   "IL1beta.Response","CD8.T.effector","C.gas_pathway","Antigen.processing.machinery"))


Fig3A<-ggplot(data = melted_cormat, aes(ord.x, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson partial correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+ 
  geom_text(aes(label=melted_p.stars), color="black", size=5)

Fig3A + coord_flip()
############################################################################################################################

#           Calculate correlation and  Plot Figure 3B 

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


all.equal(colnames(gsva_TLS_gene), rownames(mat3_transpose))
mat3_transpose_TLS<-data.frame(mat3_transpose,t(gsva_TLS_gene) )
  

summary(lm(mat3_transpose_TLS$TLS_genes~ mat3_transpose$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))
reg2_TLS_IKZF1<-ggplotRegression(lm(mat3_transpose_TLS$TLS_genes~ mat3_transpose$IKZF1+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))
reg2_TLS_IKZF1_update<-reg2_TLS_IKZF1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                        size=14),axis.text.y = element_text(face="bold",                                                                                                                                                                      size=14)) +
  xlab("Normalized expression IKZF1") + ylab("TLS signature") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


reg1_TLS_TE<-ggplotRegression(lm(mat3_transpose_TLS$TLS_genes~ mat3_transpose$TE+as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))+ as.numeric(as.factor(Ordered_pheno$batch))))


reg1_TLS_TE_update<-reg1_TLS_TE+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                 size=14),axis.text.y = element_text(face="bold",
                                                                                                                                                                                     size=14)) +
  xlab("TE") + ylab("TLS signature") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

ggarrange(reg2_TLS_IKZF1_update, reg1_TLS_TE_update,
          ncol = 1, nrow = 2, align = "v")

############################################################################################################################

#             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR TE score and IKZF1; these 2 figures make Figure 3C

############################################################################################################################

#################. Calculating KM for TE score; ; Figure 3C
#check order

all.equal(rownames(t(gsva_for_TE__genes)),Ordered_pheno$CMO_ID)

For_KM_TEscore<-data.frame(cbind(t(gsva_for_TE__genes),Ordered_pheno))

TE_features.cut <- surv_cutpoint(
  data=For_KM_TEscore,
  time = "PFS_months",
  event = "Alive.dead",
  variables = "TE_features"
)



For_KM_TEscore$TE_features_levels<-ifelse(For_KM_TEscore$TE_features> summary(TE_features.cut)[1]$cutpoint, "high TE score","low TE score")
median(For_KM_TEscore[For_KM_TEscore$TE_features_levels %in% "high TE score",]$PFS_months)
median(For_KM_TEscore[For_KM_TEscore$TE_features_levels %in% "low TE score",]$PFS_months)


sfit_TE <- survfit(Surv(PFS_months, Alive.dead)~TE_features_levels , data=For_KM_TEscore)
sfit_TE
summary(sfit_TE)

ggsurvplot(sfit_TE,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
           legend.labs=c("high TE score","low TE score"),font.legend=14,
           font.x=14,font.y=14,font.tickslab=14) 


fit_TEscore <- coxph(Surv(PFS_months, Alive.dead)~TE_features_levels + Abbreviation.for.Figures , data=For_KM_TEscore)
summary(fit_TEscore)


###ORR in TE groups
table(For_KM_TEscore$Response2,For_KM_TEscore$TE_features_levels)

fisher.test(table(For_KM_TEscore$Response2,For_KM_TEscore$TE_features_levels))

#################. Calculating KM for IKZF1; Figure 3C

all.equal(colnames(normalized_counts),Ordered_pheno$CMO_ID)

For_KM_IKZF1<-data.frame(cbind(normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),],Ordered_pheno))
colnames(For_KM_IKZF1)[1]<-"IKZF1"

IKZF1.cut <- surv_cutpoint(
  data=For_KM_IKZF1,
  time = "PFS_months",
  event = "Alive.dead",
  variables = "IKZF1"
)



For_KM_IKZF1$IKZF1_levels<-ifelse(For_KM_IKZF1$IKZF1> summary(IKZF1.cut)[1]$cutpoint, "high IKZF1","low IKZF1")

median(For_KM_IKZF1[For_KM_IKZF1$IKZF1_levels %in% "high IKZF1",]$PFS_months)
median(For_KM_IKZF1[For_KM_IKZF1$IKZF1_levels %in% "low IKZF1",]$PFS_months)

sfit_IKZF1 <- survfit(Surv(PFS_months, Alive.dead)~IKZF1_levels , data=For_KM_IKZF1)
sfit_IKZF1
summary(sfit_IKZF1)

ggsurvplot(sfit_IKZF1,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
           legend.labs=c("high IKZF1","low IKZF1"),font.legend=14,
           font.x=14,font.y=14,font.tickslab=14) 


fit_IKZF1 <- coxph(Surv(PFS_months, Alive.dead)~IKZF1_levels + Abbreviation.for.Figures, data=For_KM_IKZF1)
summary(fit_IKZF1)




###ORR in IKZF1 groups
table(For_KM_IKZF1$Response2,For_KM_IKZF1$IKZF1_levels)

fisher.test(table(For_KM_IKZF1$Response2,For_KM_IKZF1$IKZF1_levels))


##B-cells vs IKZF1 

all.equal(rownames(For_KM_IKZF1), rownames(mat_for_heatmap))


summary(lm(mat_for_heatmap$B.cell~For_KM_IKZF1$IKZF1 +as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))
reg2_Bcell_IKZF1<-ggplotRegression(lm(mat_for_heatmap$B.cell~For_KM_IKZF1$IKZF1 + as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`))))
reg2_Bcell_IKZF11_update<-reg2_Bcell_IKZF1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold",
                                                                                                                                                             size=14),axis.text.y = element_text(face="bold",                                                                                                                                                                      size=14)) +
  xlab("Normalized expression IKZF1") + ylab("B cells") +
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))



sessionInfo()

# to create a report from this script just open Rstudio and render the script using rmarkdown, note: the script and data must be all in the same folder, you first need to run Script_1.R and Script_2.R
#The command to render is rmarkdown::render("Script_3.R")


