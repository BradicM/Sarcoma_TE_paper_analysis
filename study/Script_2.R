#Author Martina Bradic
#This script requires that you first run Script_1.R, as it utilizes normalized counts data, an organized clinical data file, 
#and a normalized intergenic TE matrix. The script performs GLMnet analysis, evaluating all different models depicted in 
#Figure 2A. It also plots significant features from those models in Figure 2B.
#Additionally, it plots violin plots of normalized counts for four examples of significant 
#features that differ between immune hot and cold, as shown in Figure 2C.
#The script also performs a GLM test to assess the association between IKZF,Te score, and Immune types in the model, 
#while adjusting for batch and histology. 
#This analysis is reported in the manuscript. Furthermore, the script conducts conditional independence 
#(mutual information) tests to identify causal relationships between TEs, IKZF1,
#and the immune-hot/-cold phenotype, as reported in the manuscript.
#Calculation of correlations between TE features (Supplemental Figure 6B)  & TE score calculation are also performed by this script


library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)  
library(data.table)
library(ComplexHeatmap)
library(splineTimeR)
library(reshape2)
library(ggsignif)
library(gg.layers)
library(ggalt)
library(glmnet)
library(ggpubr)
library(bnlearn)
library(corrplot)
library(GSVA)


#FUNCTION TO NEGATE %in%
"%nin%" = Negate("%in%")



#GET function for running and extracting results for glmnet
source("GLMnet_Functions.R")


############################################################################################################################

#          RUNING AND EVALUATING FIVE DIFFERENT GLMNET MODELS 

############################################################################################################################

#GEt the full list of epigenetic genes that will be needed for the analysis
Epigentic_pathway_genes<-fread("Epigentic_pathway_genes")
Epigentic_pathway_genes<-Epigentic_pathway_genes[,-c(1)]

#539 genes are expressed
normalized_counts_all_EPI_genes<-normalized_counts[rownames(normalized_counts) %in% Epigentic_pathway_genes$Epigentic_genes,]



##################################################################################################################  
################### ####                      For cellularity model #############################################################################

set.seed(135) # this will keep the result constant, reproducible 

predictors_no_TEs<-cbind(as.numeric(as.factor(Ordered_pheno$`Abbreviation for Figures`)), as.numeric(as.factor(Ordered_pheno$batch)))
colnames(predictors_no_TEs)<-c("Histology","batch")
rownames(predictors_no_TEs)<-Ordered_pheno$CMO_ID

predictors_no_TEs<-data.frame(predictors_no_TEs)

Cellularity_model<-glmnet_binary_function1000x(as.matrix(predictors_no_TEs),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))
Cellularity_model<-data.frame(Cellularity_model[[1]])
Cellularity_model_features<-Function_for_glmnet(as.matrix(predictors_no_TEs),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))
Cellularity_model_features$glmnet_results_df

#Plot features that contribute to the model 
features<-data.frame(Cellularity_model_features$glmnet_results_df)
features$names<-rownames(features)

#############################################################################

################### ####      For TE model  #############################################################################

#############################################################################
set.seed(135)
TE_matrix<-INTERGENIC_TE_normalized_expression_1052_repeats
TE_matrix<-t(TE_matrix)
#Check if the ORder pheno is the same sample order with TE_matrix
all.equal(rownames(TE_matrix) , Ordered_pheno$CMO_ID)


predictors<-data.frame(predictors_no_TEs,TE_matrix)


TE_and_cellularity_model<-glmnet_binary_function1000x(as.matrix(predictors),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))
TE_and_cellularity_model<-data.frame(TE_and_cellularity_model[[1]])


TE_and_Cellularity_model_features<-Function_for_glmnet(as.matrix(predictors),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))

TE_and_Cellularity_model_features$glmnet_results_df

features2<-data.frame(TE_and_Cellularity_model_features$glmnet_results_df)
features2$names<-rownames(features2)
features2<-features2[-c(1),]
ggplot(features2, 
       aes(x=names, y=s1, label = round(s1, 1))) + 
  geom_lollipop(point.size = 3, point.colour = "cadetblue") +
  #geom_text(nudge_y = 5) +
  coord_flip() +
  #theme_minimal() +
  labs(y="Contribution", x="Feature") +
  theme(text = element_text(size=20)) 



################### ####      For TE model  shuffled #############################################################################

#############################################################################
set.seed(135) # this will keep the result constant, reproducable 

TE_and_cellularity_model_shuffeled<-glmnet_shuffeled_function1000x_get_all_1000R2(as.matrix(predictors),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))
TE_and_Cellularity_model_features_shuffeled<-Function_for_glmnet(as.matrix(predictors),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))

TE_and_Cellularity_model_features_shuffeled$glmnet_results_df


#############################################################################

################### ####      For  celularity + all EPI GENES  #############################################################################

#############################################################################
set.seed(135) # this will keep the result constant, reproducable 

predictors_EPI<-data.frame(predictors_no_TEs,t(normalized_counts_all_EPI_genes))

#Check if the ORder pheno is the same sample order with TE_matrix
all.equal(rownames(TE_matrix) , Ordered_pheno$CMO_ID)
all.equal(rownames(predictors_EPI),Ordered_pheno$CMO_ID)


Cellularity_and_epig_gene<-glmnet_binary_function1000x(as.matrix(predictors_EPI),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))
Cellularity_and_epig_gene_features<-Function_for_glmnet(as.matrix(predictors_EPI),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))


#Plot features that contribute to the model 
features5<-data.frame(Cellularity_and_epig_gene_features$glmnet_results_df)

features5$names<-rownames(features5)
features5<-features5[-c(1),]
ggplot(features5, 
       aes(x=names, y=s1, label = round(s1, 1))) + 
  geom_lollipop(point.size = 3, point.colour = "cadetblue") +
  #geom_text(nudge_y = 5) +
  coord_flip() +
  #theme_minimal() +
  labs(y="Contribution", x="Feature") +
  theme(text = element_text(size=20)) 


################### ####      For  cellularity + all EPigenetic genes  shuffled #############################################################################
Cellularity_and_epig_gene_shuffeled<-glmnet_shuffeled_function1000x_get_all_1000R2(as.matrix(predictors_EPI),as.matrix(as.numeric(factor(Ordered_pheno$Immune_subtypes))))


############################################################################################################################

#           Plotting Fig 2A and B, feature from TE model and Epi model 

############################################################################################################################

#Plotting differnet model contribution Figure 2A
Results_models_EPi<-data.frame(Cellularity_model, TE_and_cellularity_model,t(TE_and_cellularity_model_shuffeled), Cellularity_and_epig_gene,t(Cellularity_and_epig_gene_shuffeled))
colnames(Results_models_EPi)<-c("Basic model", "Basic model + TE ", "Basic model + TE schuffeled","Basic model + EPI genes","Basic model + EPI genes schuffeled")


Results_models_EPi_2<-melt(Results_models_EPi)
colnames(Results_models_EPi_2)<-c("model","r2")


# function for computing mean, DS, max and min values


min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
  
}

p2 <- ggplot(aes(y = r2, x = factor(model)), data = Results_models_EPi_2)
p2 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot2") +  
  xlab("Model") + ylab("R2") +
  theme_bw()+
  theme(axis.line.x = element_blank(), 
        axis.line.y = element_blank(),
        axis.text=element_text(size=5),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(angle=45,hjust=1)) +
  #geom_text(aes(label = variable), size=7,
  #        position = position_stack(vjust = 0.5)) +
  theme(axis.text = element_text(colour = "black", size = rel(1.3)),panel.grid.major = element_blank(), panel.grid.minor = element_blank())


####Test statistics for shuffled vs unshuffled model 
 t.test(Results_models_EPi_2[Results_models_EPi_2$model %in% "Basic model + TE ",]$r2, Results_models_EPi_2[Results_models_EPi_2$model %in% "Basic model + TE schuffeled",]$r2)
 t.test(Results_models_EPi_2[Results_models_EPi_2$model %in% "Basic model + EPI genes",]$r2, Results_models_EPi_2[Results_models_EPi_2$model %in% "Basic model + EPI genes schuffeled",]$r2)
 
 
 #plotting features from significant models Figure 2B
 
 features_2_5<-data.frame(rbind(features2,features5))
 ggplot(features_2_5, 
        aes(x=names, y=s1, label = round(s1, 1))) + 
   geom_lollipop(point.size = 3, point.colour = "cadetblue") +
   #geom_text(nudge_y = 5) +
   coord_flip() +
   #theme_minimal() +
   labs(y="Contribution", x="Feature") +
   theme(text = element_text(size=20)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
 
 
 
 
 ############################################################################################################################
 
 #           Plotting Fig 2C
 
 ############################################################################################################################
 
 #violin plots of normalized counts for 4 examples of significant features that differ between immune hot and cold, Figure 2C
 
 my_colors= c('blue','red' )
 
 
 all.equal(colnames(normalized_counts),colnames(INTERGENIC_TE_normalized_expression_1052_repeats))
 mat1<-rbind(normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),],INTERGENIC_TE_normalized_expression_1052_repeats[rownames(INTERGENIC_TE_normalized_expression_1052_repeats) %in% c("HERVL74.int","MER45A", "MER61F","MER57F","Tigger17a"),])
 mat1<-data.frame(t(mat1))
 rownames(mat1)<-gsub (".","-",fixed=TRUE,rownames(mat1))
 all.equal(rownames(mat1),Ordered_pheno$CMO_ID)
 colnames(mat1)[1]<-"IKZF1"
 
 p1A<-ggplot(mat1, aes(x=Ordered_pheno$Immune_subtypes, y=IKZF1,fill = Ordered_pheno$Immune_subtypes)) + 
   geom_violin(alpha = 0.5) +
    #ylim(0,NA) +
   scale_fill_manual(values = my_colors) +
   geom_point(position = position_jitter(seed = 1, width = 0.2)) +
   theme(legend.position = "none") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),
         axis.title.x = element_blank(),axis.title.y = element_text(size=12), 
         axis.text.x = element_text( size=10),axis.text.y = element_text(size=12),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   geom_signif(comparisons = list(c("Immune cold", "Immune hot")), 
               map_signif_level=TRUE) +
  ggtitle("IKZF1") + ylab("Normalized expression") 
   
  
 p2A<-ggplot(mat1, aes(x=Ordered_pheno$Immune_subtypes, y=MER57F,fill = Ordered_pheno$Immune_subtypes)) + 
   geom_violin(alpha = 0.5) +
   # ylim(0,NA) +
   scale_fill_manual(values = my_colors) +
   geom_point(position = position_jitter(seed = 1, width = 0.2)) +
   theme(legend.position = "none") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),
         axis.title.x = element_blank(),axis.title.y = element_text(size=12), 
         axis.text.x = element_text( size=10),axis.text.y = element_text(size=12),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   geom_signif(comparisons = list(c("Immune cold", "Immune hot")), 
               map_signif_level=TRUE) +
   ggtitle("MER57F") + ylab("Normalized expression") 
 
 p3A<-ggplot(mat1, aes(x=Ordered_pheno$Immune_subtypes, y=MER61F,fill = Ordered_pheno$Immune_subtypes)) + 
   geom_violin(alpha = 0.5) +
  #  ylim(0,NA) +
   scale_fill_manual(values = my_colors) +
   geom_point(position = position_jitter(seed = 1, width = 0.2)) +
   theme(legend.position = "none") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),
         axis.title.x = element_blank(),axis.title.y = element_text(size=12), 
         axis.text.x = element_text( size=10),axis.text.y = element_text(size=12),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   geom_signif(comparisons = list(c("Immune cold", "Immune hot")), 
               map_signif_level=TRUE) +
   ggtitle("MER61F") + ylab("Normalized expression") 
 
 p4A<-ggplot(mat1, aes(x=Ordered_pheno$Immune_subtypes, y=MER45A,fill = Ordered_pheno$Immune_subtypes)) + 
   geom_violin(alpha = 0.5) +
   # ylim(0,NA) +
   scale_fill_manual(values = my_colors) +
   geom_point(position = position_jitter(seed = 1, width = 0.2)) +
   theme(legend.position = "none") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),
         axis.title.x = element_blank(),axis.title.y = element_text(size=12), 
         axis.text.x = element_text( size=10),axis.text.y = element_text(size=12),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   geom_signif(comparisons = list(c("Immune cold", "Immune hot")), 
               map_signif_level=TRUE) +
   ggtitle("MER45A") + ylab("Normalized expression") 
 
 
 
 ggarrange(p1A, p2A,p3A,p4A,
           ncol = 2, nrow = 2, align = "v")
 
 
####################Calculate correlations between TE features (Supplemental Figure 6B)  & TE score calculation#######################################
 
 TE_features<-c(rownames(TE_and_Cellularity_model_features$glmnet_results_df)[-c(1)])
 TE_features<-data.frame(TE_features)
 
 
 
 TE_genes_df<-as.data.frame(TE_matrix[,colnames(TE_matrix) %in% as.character(TE_features$TE_features)])
 M_TE_corr <-cor(TE_genes_df)
 
 
 
 corrplot(M_TE_corr, type="upper",method = 'circle', order = 'AOE',tl.col = "black",tl.cex = 0.8,
          col=colorRampPalette(c("blue","white","red"))(200))
 
 
 #calculate TE scores
 geneSets<-as.list(TE_features)
 geneSets$TE_features<-as.character(geneSets$TE_features)
 geneSets<-lapply(geneSets, function(z){ z[!is.na(z) & z != ""]})
 
 TE_genes_matrix<-as.matrix(TE_matrix[,colnames(TE_matrix) %in% as.character(TE_features$TE_features)])
 
 gsva_for_TE__genes <-gsva(t(TE_genes_matrix),geneSets,method="zscore", kcdf="Gaussian",verbose=FALSE)
 
 
 ########################### glm to test association between IKZF, Te score and Immune types and conditional independence #######
 all.equal(Ordered_pheno$CMO_ID,colnames(gsva_for_TE__genes),colnames(normalized_counts_all_EPI_genes))
 
 
 #Linear model to test association between TE score, and IKZF1 expression with immune type 
 
 table_for_model<-data.frame(Ordered_pheno$Immune_subtypes , gsva_for_TE__genes["TE_features",] ,normalized_counts_all_EPI_genes["IKZF1",], Ordered_pheno$`Abbreviation for Figures` , Ordered_pheno$batch)
 colnames(table_for_model)<-c("Immune_subtypes","TE.score","IKZF1","Abbreviation for Figures","batch")
 
 fit<-glm(as.numeric(as.factor(Immune_subtypes)) ~  TE.score + IKZF1 + `Abbreviation for Figures` + batch, data=table_for_model)
 
 summary(fit)
 
 
 #We need to convert immune type into numerical value for bnlearn to work 

 #hypotheses 1
 hypo1<-ci.test(y=as.numeric(as.factor(table_for_model$Immune_subtypes)), x=table_for_model[,"TE.score"], z=table_for_model[,"IKZF1"])
 
 #hypotheses 2
 hypo2<-ci.test(y=as.numeric(as.factor(table_for_model$Immune_subtypes)), x=table_for_model[,"IKZF1"], z= table_for_model[,"TE.score"])
 
 #hypotheses 3
 hypo3<-ci.test(x=table_for_model[,"TE.score"], y=table_for_model[,"IKZF1"],z=as.numeric(as.factor(table_for_model$Immune_subtypes)))
 
 
 all_3_poss<-data.frame(hypo1$p.value,hypo2$p.value,hypo3$p.value)
 print(all_3_poss)
 
 
 
 sessionInfo()
 
 # to create a report from this script just open Rstudio and render the script using rmarkdown, note: the script and data must be all in the same folder, you first need to run Script_1.R
 #The command to render is rmarkdown::render("Script_2.R")
 
