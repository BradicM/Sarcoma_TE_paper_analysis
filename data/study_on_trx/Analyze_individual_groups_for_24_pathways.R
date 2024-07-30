



get_individual_pairs_comparison_calculations_24_pathways<- function(Pairs_comparison){
  
  # choose the cohort
  
  Ordered_pheno_preselection<-Ordered_pheno[Ordered_pheno$CMO_ID2 %in% Pairs_comparison$CMOID,]
  
  {
    
    ##########################################. RUN PAIRS COMPARISON #########################################.
    
    #select pairs for the selection above
    identified_pairs<-table(Ordered_pheno_preselection$CMO_ID2) ==2
    pair_names<-names(identified_pairs[identified_pairs==TRUE])
    
    #this is the final table that will be used for the comparison, this has only paired samples selected from above criteria; Ordered_pheno_preselection
    Ordered_pheno_selection<-Ordered_pheno_preselection[Ordered_pheno_preselection$CMO_ID2 %in% pair_names,]
    
    
    
    
    #########################################.  Calculate difference between baseline and on-trx for each cell type 
    
   
    
    
    gsva_24_gene_lists_selection<-gsva_24_gene_lists_df[gsva_24_gene_lists_df$CMO_ID %in% Ordered_pheno_selection$CMO_ID,]
    all.equal(gsva_24_gene_lists_selection$CMO_ID, Ordered_pheno_selection$CMO_ID)
    
    gsva_24_gene_ordered_selection_with_pheno<-left_join(gsva_24_gene_lists_selection,Ordered_pheno_selection)
    head(gsva_24_gene_ordered_selection_with_pheno[,1:25])
    
    result_t.test_immune_cells = matrix(NA,nrow = 25)
    Immune_cell<-colnames(gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline",])[1:25]
    
    for(i in 1:length(Immune_cell)) {
      test_results<-t.test(gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline",][,i],gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx",][,i], paired = TRUE)
      
      result_t.test_immune_cells[i,] <- test_results$p.value
    }
    result_t.test_immune_cells<-setNames(data.frame(result_t.test_immune_cells), c("p_values"))
    rownames(result_t.test_immune_cells)<-Immune_cell
    # print(result_t.test_immune_cells)
    
    #t.test(gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline" & gsva_24_gene_ordered_selection_with_pheno$Response2 %in% "Responder",][,Immune_cell],gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx" & gsva_24_gene_ordered_selection_with_pheno$Response2 %in% "Responder",][,Immune_cell],paired=TRUE)
    #t.test(gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline" & gsva_24_gene_ordered_selection_with_pheno$Response2 %in% "Nonreseponder",][,Immune_cell],gsva_24_gene_ordered_selection_with_pheno[gsva_24_gene_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx" & gsva_24_gene_ordered_selection_with_pheno$Response2 %in% "Nonreseponder",][,Immune_cell],paired=TRUE)
    
    
    
    plot_interferon<-ggplot(gsva_24_gene_ordered_selection_with_pheno, aes(Sample_Timepoint, Type.II.IFN.Response, fill=Sample_Timepoint)) + 
      geom_violin()+ 
      
      # geom_point() is used to make points at data values 
      geom_point()+ 
      # geom_line() joins the paired datapoints 
      geom_line(aes(group=CMO_ID2,color = Response2))  +
      scale_color_manual(values=c("blue", "red")) +
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14)) +
      
      scale_fill_manual(values = c(Baseline = 'gray', 
                                   Ontrx = 'white'))
    
    
    plot_CD8<-ggplot(gsva_24_gene_ordered_selection_with_pheno, aes(Sample_Timepoint, CD8.T.effector, fill=Sample_Timepoint)) + 
      geom_violin()+ 
      
      # geom_point() is used to make points at data values 
      geom_point()+ 
      # geom_line() joins the paired datapoints 
      geom_line(aes(group=CMO_ID2,color = Response2))  +
      scale_color_manual(values=c("blue", "red")) +
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14)) +
      
      scale_fill_manual(values = c(Baseline = 'gray', 
                                   Ontrx = 'white'))
    
    
    plot_NFkB<-ggplot(gsva_24_gene_ordered_selection_with_pheno, aes(Sample_Timepoint, NFkB.response, fill=Sample_Timepoint)) + 
      geom_violin()+ 
      
      # geom_point() is used to make points at data values 
      geom_point()+ 
      # geom_line() joins the paired datapoints 
      geom_line(aes(group=CMO_ID2,color = Response2))  +
      scale_color_manual(values=c("blue", "red")) +
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14)) +
      
      scale_fill_manual(values = c(Baseline = 'gray', 
                                   Ontrx = 'white'))
    
    
    plot_C_GAS<-ggplot(gsva_24_gene_ordered_selection_with_pheno, aes(Sample_Timepoint, C_GAS_gene, fill=Sample_Timepoint)) + 
      geom_violin()+ 
      
      # geom_point() is used to make points at data values 
      geom_point()+ 
      # geom_line() joins the paired datapoints 
      geom_line(aes(group=CMO_ID2,color = Response2))  +
      scale_color_manual(values=c("blue", "red")) +
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14)) +
      
      scale_fill_manual(values = c(Baseline = 'gray', 
                                   Ontrx = 'white'))
   
    
  }

  list(result_t.test_immune_cells,plot_interferon,plot_CD8,plot_NFkB,plot_C_GAS)
  
}
