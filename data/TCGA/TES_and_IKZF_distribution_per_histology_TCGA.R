all.equal(colnames(normalized_counts),colnames(INTERGENIC_TE_normalized_expression_1052_repeats))
TE_intergenic_for_heatmap_1<-rbind(normalized_counts[rownames(normalized_counts) %in% c("IKZF1"),],INTERGENIC_TE_normalized_expression_1052_repeats[rownames(INTERGENIC_TE_normalized_expression_1052_repeats) %in% c("MER57F","MER61F","MER44A","LTR29","AluYk3","FordPrefect_a"),])

TE_intergenic_for_heatmap_df<-data.frame(t(TE_intergenic_for_heatmap_1))
TE_intergenic_for_heatmap_df$patient<-rownames(TE_intergenic_for_heatmap_df)
TE_intergenic_for_heatmap_df<-left_join(TE_intergenic_for_heatmap_df,Ordered_pheno)
TE_intergenic_for_heatmap_df_sub<-TE_intergenic_for_heatmap_df[TE_intergenic_for_heatmap_df$patient %in% TE_intergenic_for_heatmap_df[TE_intergenic_for_heatmap_df$`short histo` %in% c("LMS","DDLPS","UPS"),]$patient,]
colnames(TE_intergenic_for_heatmap_df_sub)[1]<-"IKZF1"
TE_intergenic_for_heatmap_df_sub2<-TE_intergenic_for_heatmap_df_sub[,colnames(TE_intergenic_for_heatmap_df_sub) %in% c("MER57F","MER61F","MER44A","LTR29","AluYk3","FordPrefect_a","IKZF1","short histo","purity")]


## AluYk3          0.230245986
## FordPrefect_a   0.183558306
## LTR103b_Mam     0.029486498
## LTR16E2         0.114007700
## LTR1C1          0.049637768
## LTR29           0.240670883
## LTR33           0.032326698
## LTR9            0.046238168
## MER101         -0.007169067
## MER44A          0.219633580
## MER57F          0.174049127
## MER61F          0.178620263
## hAT.N1_Mam      0.033856884

#MER57F,MER61F,MER44A,LTR29,AluYk3,FordPrefect_a

sig_TEs<-c("MER45A","MER57F","Tigger17a","MER61F","LTR104_Mam","HERVL74.int")
my_colors = c("#DD67C0","#F4A460", "#0000FF","#82D992")
IKZF1<-ggplot(TE_intergenic_for_heatmap_df_sub2, aes(`short histo`, IKZF1), fill=`Abbreviation for Figures`) + 
  
geom_violin(alpha = 0.5) +
  #ylim(0,NA) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),
        axis.title.x = element_blank(),axis.title.y = element_text(size=12), 
        axis.text.x = element_text( size=10),axis.text.y = element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
geom_signif(comparisons = list(c("DDLPS","UPS"),c("LMS","UPS"),c("DDLPS","UPS")), 
 #           geom_signif(comparisons = list(c("ANGS","UPS"),c("LMS","UPS"),c("LPS","UPS")), 
                        
                      map_signif_level=TRUE)
 #


ggarrange(TE1, TE2,TE3,TE4,TE5,TE6,
          ncol = 2, nrow = 3, align = "v")
