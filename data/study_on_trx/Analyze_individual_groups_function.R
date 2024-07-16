get_individual_pairs_comparison_calculations<- function(Pairs_comparison){
  
  # choose the cohort
  
  Ordered_pheno_preselection<-Ordered_pheno[Ordered_pheno$CMO_ID2 %in% Pairs_comparison$CMOID,]
  
  {
    
    ##########################################. RUN PAIRS COMPARISON #########################################.
    
    #select pairs for the selection above
    identified_pairs<-table(Ordered_pheno_preselection$CMO_ID2) ==2
    pair_names<-names(identified_pairs[identified_pairs==TRUE])
    
    #this is the final table that will be used for the comparison, this has only paired samples selected from above criteria; Ordered_pheno_preselection
    Ordered_pheno_selection<-Ordered_pheno_preselection[Ordered_pheno_preselection$CMO_ID2 %in% pair_names,]
    
    
    
    
    ##################################################.  Anlsysis , Firt plot heatmap of immune proportion of selected group and genes 
    
    
    mat_for_heatmap_ordered_selection<-t(mat_immuno2_no_batch[,colnames(mat_immuno2_no_batch) %in% Ordered_pheno_selection$CMO_ID])
    
    #order pheno and matrix
    mat_for_heatmap_ordered_selection<-mat_for_heatmap_ordered_selection[order(rownames(mat_for_heatmap_ordered_selection)),]
    Ordered_pheno_selection<-Ordered_pheno_selection[order(Ordered_pheno_selection$CMO_ID),]
    
    ann_for_selection <- data.frame(Ordered_pheno_selection[,c("Response","Abbreviation for Figures","Immune_subtypes","Sample_Timepoint")]) 
    rownames(ann_for_selection)<-Ordered_pheno_selection$CMO_ID
    colnames(ann_for_selection) <- c("Response","Subtype","Immune_subtypes","Sample_Timepoint") 
    
    colours_selection <- list(
      'Response' = c("PD"="#543005","SD"="#A6611A" , "CR/PR"= "#F6E8C3"),
      'Subtype'=c ("ANGS"="#DD67C0" ,"CHS"="#A848E0","DDLS"= "#DEABCA","LMS"= "#F4A460","OS"= "#8881D2" ,"SARCNOS"="#8DB6D1","UPS"= "#82D992","EHE"="#FFFF00","LPS"="#0000FF","MFS"="#B22222","Other"="#FF0000","SBRC"="#85E0D6","ASPS"="grey"),
      'Immune_subtypes'=c("Immune cold"="blue","Immune hot"="red"),
      'Sample_Timepoint'=c("Baseline"="darkgreen","Ontrx"="darkred")
    )
    
    
    
    colAnn_selection <- HeatmapAnnotation(df = ann_for_selection,
                                          which = "col",
                                          col = colours_selection,
                                          annotation_width = unit(c(1, 4), "cm"),
                                          gap = unit(1,"mm"))
    
    
    all.equal(rownames(mat_for_heatmap_ordered_selection),Ordered_pheno_selection$CMO_ID)
    all.equal(rownames(mat_for_heatmap_ordered_selection),rownames(ann_for_selection))
    
    Heat1_A_B_C_selection<-Heatmap(t(mat_for_heatmap_ordered_selection), name = "Z-score",row_gap = unit(5, "mm"),clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_selection),column_split =Ordered_pheno_selection$CMO_ID2,
                                   top_annotation =colAnn_selection,show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE)
    
    
    
    IKZF1<-data.frame(t(normalized_counts[rownames(normalized_counts) %in% c("IKZF1", "CTLA4","CD274","LAG3"),]))
    
    
    mat_for_heatmap_ordered_selection_with_pheno<-data.frame(mat_for_heatmap_ordered_selection)
    mat_for_heatmap_ordered_selection_with_pheno$CMO_ID<-rownames(mat_for_heatmap_ordered_selection)
    mat_for_heatmap_ordered_selection_with_pheno<-left_join(mat_for_heatmap_ordered_selection_with_pheno,Ordered_pheno_selection)
    
    IKZF1_selection<-IKZF1[rownames(IKZF1) %in% mat_for_heatmap_ordered_selection_with_pheno$CMO_ID,]
    
    #order IKZF1
    IKZF1_selection<-IKZF1_selection[order(rownames(IKZF1_selection)),]
    
    
    all.equal(rownames(IKZF1_selection),Ordered_pheno_selection$CMO_ID)
    all.equal(rownames(IKZF1_selection),mat_for_heatmap_ordered_selection_with_pheno$CMO_ID)
    
    
    
    
    #IKZF1_selection_no_batch<-removeBatchEffect(IKZF1_selection, IKZF1_selection$batch,as.numeric(asIKZF1_selection$`Abbreviation for Figures`))
    
    ht1C_only_responder<-Heatmap(t(IKZF1_selection[,c("CD274","CTLA4","IKZF1","LAG3")]),  name = "ICI genes",clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D" ,row_names_gp = gpar(fontsize = 10),column_labels = rownames(ann_for_selection),column_split =Ordered_pheno_selection$CMO_ID2,
                                 show_column_dend = FALSE,show_row_dend = FALSE,show_column_names = FALSE)
    
    
    FIG_selection_combined = Heat1_A_B_C_selection %v% ht1C_only_responder 
    
    
    
    #########################################.  Calculate difference between baseline and on-trx for each cell type 
    
    mat_for_heatmap_ordered_selection_with_pheno<-cbind(mat_for_heatmap_ordered_selection,Ordered_pheno_selection)
    
    
    result_t.test_immune_cells = matrix(NA,nrow = 11)
    Immune_cell<-colnames(mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline",])[1:11]
    
    for(i in 1:length(Immune_cell)) {
      test_results<-t.test(mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline",][,i],mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx",][,i], paired = TRUE)
      
      result_t.test_immune_cells[i,] <- test_results$p.value
    }
    result_t.test_immune_cells<-setNames(data.frame(result_t.test_immune_cells), c("p_values"))
    rownames(result_t.test_immune_cells)<-Immune_cell
   # print(result_t.test_immune_cells)
    
    #t.test(mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline" & mat_for_heatmap_ordered_selection_with_pheno$Response2 %in% "Responder",][,Immune_cell],mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx" & mat_for_heatmap_ordered_selection_with_pheno$Response2 %in% "Responder",][,Immune_cell],paired=TRUE)
    #t.test(mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Baseline" & mat_for_heatmap_ordered_selection_with_pheno$Response2 %in% "Nonreseponder",][,Immune_cell],mat_for_heatmap_ordered_selection_with_pheno[mat_for_heatmap_ordered_selection_with_pheno$Sample_Timepoint %in% "Ontrx" & mat_for_heatmap_ordered_selection_with_pheno$Response2 %in% "Nonreseponder",][,Immune_cell],paired=TRUE)
    
    
    
    
    ####################################################.   calculate IKZF1 changes between baseline and ontrx
    IKZF1_selection<-cbind(IKZF1_selection,Ordered_pheno_selection)
    plotIKZF1_select<-ggplot(IKZF1_selection, aes(Sample_Timepoint, IKZF1, fill=Sample_Timepoint)) + 
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
    
    
    
    result_t.test_IKZF1<-t.test(IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Baseline",]$IKZF1,IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Ontrx",]$IKZF1,paired=TRUE)
    #t.test(IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Baseline" & IKZF1_selection$Response2 %in% "Responder",]$IKZF1,IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Ontrx" & IKZF1_selection$Response2 %in% "Responder",]$IKZF1,paired=TRUE)
    #t.test(IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Baseline" & IKZF1_selection$Response2 %in% "Nonreseponder",]$IKZF1,IKZF1_selection[IKZF1_selection$Sample_Timepoint %in% "Ontrx" & IKZF1_selection$Response2 %in% "Nonreseponder",]$IKZF1,paired=TRUE)
    
    
    
     
    
    ##########################################################################################  calculate TEscore changes between baseline and ontrx 
    
    TE_genes_matrix_df_selection<-TE_genes_matrix_df[TE_genes_matrix_df$CMO_ID2 %in% Ordered_pheno_selection$CMO_ID2,]
    
    
    plotTEscore_select<-ggplot(TE_genes_matrix_df_selection, aes(Sample_Timepoint, TE_features, fill=Sample_Timepoint)) + 
      
      geom_violin()+ 
      
      # geom_point() is used to make points at data values 
      geom_point() +
      
      # geom_line() joins the paired datapoints 
      geom_line(aes(group=CMO_ID2,color = Response2))  +
      scale_color_manual(values=c("blue", "red")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", size=14),axis.text.y = element_text(face="bold", size=14)) +
      
      scale_fill_manual(values = c(Baseline = 'gray', 
                                   Ontrx = 'white'))
    
    
    
    
    #for all
    result_t.test_TE<- t.test(TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Baseline",]$TE_features,TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Ontrx",]$TE_features,paired=TRUE)
    
    #for responder
    #t.test(TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Baseline" & TE_genes_matrix_df_selection$Response2 %in% "Responder",]$TE_features,TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Ontrx" & TE_genes_matrix_df_selection$Response2 %in% "Responder",]$TE_features,paired=TRUE)
    
    #for non-responder
    #t.test(TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Baseline" & TE_genes_matrix_df_selection$Response2 %in% "Nonreseponder",]$TE_features,TE_genes_matrix_df_selection[TE_genes_matrix_df_selection$Sample_Timepoint %in% "Ontrx" & TE_genes_matrix_df_selection$Response2 %in% "Nonreseponder",]$TE_features,paired=TRUE)
    
    
    
    ###join 2 plots together 
    
    # make this for each selection 
    
    boxplots_select<-ggarrange(plotIKZF1_select, plotTEscore_select,ncol = 1, nrow = 2, align = "v")
   
    
  }

  list(FIG_selection_combined, result_t.test_immune_cells,boxplots_select,result_t.test_IKZF1, result_t.test_TE)
  
}
