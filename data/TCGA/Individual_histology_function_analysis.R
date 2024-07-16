#calculate statistic for each cohort separately 
# This analysis requires for Script_1_TCGA.R, Script_2_TCGA.R, and Script_3_TCGA.R to be executed


#####################################################################################################################################################
############################################   FUNCTION TO RUN Kaplan Meier analysis for individual features immune type IKZF1, TE score, TLS#######
############################################        FOR EACH HISTOLOGY ##############################################################################
####################################################################################################################################

get_individual_TCGA_histology_KM_calcualtions<- function(cohort){
  
  # choose the cohort
  
  Ordered_pheno_select<-Ordered_pheno[Ordered_pheno$`short histo` %in% cohort,]
  
  {
    
    
    ############################################################################################################################
    
    #             KAPLAN MEIEIR ANALYSIS AND PLOTS FOR Immune hot/cold for individual cphorts
    
    ############################################################################################################################
    
    
    #Model for immune hot/cold, Figure 4D
    sfit_select <- survminer::surv_fit(Surv(OS_months, OS_status)~Immune_subtypes , data=Ordered_pheno_select)
    
    plot_Immune_type<-ggsurvplot(sfit_select,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("#0000FF","red"),
                                 legend.labs=c("Immune-cold","Immune-hot"),font.legend=14,
                                 font.x=14,font.y=14,font.tickslab=14) 
    
    
    
    fit_select <- coxph(Surv(OS_months, OS_status)~Immune_subtypes  + `pathologic tumor size`, data=Ordered_pheno_select)
    Coxph_Immune_type_select<-summary(fit_select)
    
    
    
    median(Ordered_pheno_select[Ordered_pheno_select$Immune_subtypes %in% "Immune hot",]$OS_months)
    median(Ordered_pheno_select[Ordered_pheno_select$Immune_subtypes %in% "Immune cold",]$OS_months)
    
    
    
    #Model for IKZF1 high/low, figure 4F
    
    normalized_counts_select<-normalized_counts[,colnames(normalized_counts) %in% Ordered_pheno_select$patient]
    all.equal(colnames(normalized_counts_select),Ordered_pheno_select$patient)
    For_KLM_IKZF1_select<-data.frame(Ordered_pheno_select,normalized_counts_select[rownames(normalized_counts_select) %in% "IKZF1",])
    colnames(For_KLM_IKZF1_select)[70]<-"IKZF1"
    #Let's derive a cutpoint for IKZF1 expression 
    
    IKZF21.surv_rnaseq.cut_select <- surv_cutpoint(
      data=For_KLM_IKZF1_select,
      time = "OS_months",
      event = "OS_status",
      variables = "IKZF1"
    )
    summary(IKZF21.surv_rnaseq.cut_select)
    
    
    For_KLM_IKZF1_select$IKZF_high_low <-ifelse (For_KLM_IKZF1_select$IKZF1  > summary(IKZF21.surv_rnaseq.cut_select)[1]$cutpoint, "high expression","low expression")
    
    sfit_IKZF1_select <- survminer::surv_fit(Surv(OS_months, OS_status)~IKZF_high_low, data=For_KLM_IKZF1_select)
    
    plot_IKZF1<-ggsurvplot(sfit_IKZF1_select,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
                           legend.labs=c("high IKZF1","low IKZF1"),font.legend=14,
                           font.x=14,font.y=14,font.tickslab=14) 
    
    
    
    
    
    fit_IKZF1_select<- coxph(Surv(OS_months, OS_status)~IKZF1  + pathologic.tumor.size, data=For_KLM_IKZF1_select)
    Coxph_IKZF1_select<-summary(fit_IKZF1_select)
    
    
    #################. Calculating KM for TE score; 
    #check order
    gsva_for_TE__genes_select<-data.frame(t(gsva_for_TE__genes[,colnames(gsva_for_TE__genes) %in% Ordered_pheno_select$patient]))
    rownames(gsva_for_TE__genes_select)<-"TE_features"
    colnames(gsva_for_TE__genes_select)<-gsub(".","-",fixed=TRUE, colnames(gsva_for_TE__genes_select))
    all.equal(rownames(t(gsva_for_TE__genes_select)),Ordered_pheno_select$patient)
    
    For_KM_TEscore_select<-data.frame(cbind(t(gsva_for_TE__genes_select),Ordered_pheno_select))
    
    TE_features.cut_select <- surv_cutpoint(
      data=For_KM_TEscore_select,
      time = "OS_months",
      event = "OS_status",
      variables = "TE_features"
    )
    
    
    
    For_KM_TEscore_select$TE_features_levels<-ifelse(For_KM_TEscore_select$TE_features> summary(TE_features.cut_select)[1]$cutpoint, "high TE score","low TE score")
    
    sfit_TE_select <- survminer::surv_fit(Surv(OS_months, OS_status)~TE_features_levels , data=For_KM_TEscore_select)
    
    
    plot_TEscore<-ggsurvplot(sfit_TE_select,pval=TRUE, conf.int=FALSE, risk.table=TRUE,risk.table.height=.15,palette = c("red","#0000FF"),
                             legend.labs=c("high TE score","low TE score"),font.legend=14,
                             font.x=14,font.y=14,font.tickslab=14) 
    
    
    fit_TEscore_select <- coxph(Surv(OS_months, OS_status)~TE_features_levels  + pathologic.tumor.size, data=For_KM_TEscore)
    Coxph_TEscore_select<-summary(fit_TEscore_select)
    
    
    
    
  }
  #return(plot_Immune_type)
  #return(plot_IKZF1)
  #return(plot_TEscore)
  #return(Coxph_Immune_type)
  list(plot_Immune_type, Coxph_Immune_type_select,plot_IKZF1,Coxph_IKZF1_select,plot_TEscore,Coxph_TEscore_select)
  
}

