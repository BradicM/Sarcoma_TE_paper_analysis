
#Organized data

Xcell_Results_for_genes<-fread("Xcell_Result")

Xcell_Results_for_genes<-data.frame(Xcell_Results_for_genes)
rownames(Xcell_Results_for_genes)<-Xcell_Results_for_genes$V1
Xcell_Results_for_genes<-Xcell_Results_for_genes[,-c(1)]

colnames(Xcell_Results_for_genes)<-gsub(".","-",colnames(Xcell_Results_for_genes), fixed=TRUE)
colnames(Xcell_Results_for_genes)<-gsub("X","",colnames(Xcell_Results_for_genes), fixed=TRUE)


#Calculate lymphoid and myeloid cell content

lymphoid_content_cells<-c("CD8+ T-cells","NK cells","CD4+ naive T-cells","B-cells","CD4+ T-cells","CD8+ Tem","Tregs","Plasma cells","CD4+ Tcm","CD4+ Tem","Memory B-cells","CD8+ Tcm","naive B-cells","CD4+ memory T-cells","pro B-cells","Class-switched memory B-cells","Th2 cells","Th1 cells","CD8+ naive T-cells","NKT","Tgd cells")
lymphoid_content_cells_sum<-data.frame(colSums( Xcell_Results_for_genes[which(rownames(Xcell_Results_for_genes) %in% lymphoid_content_cells),] ))
myeloid_content_cells<-c("Monocytes", "Macrophages", "DC", "Neutrophils", "Eosinophils", "Macrophages M1", "Macrophages M2", "aDC", "Basophils", "cDC" , "pDC", "iDC","Mast cells")
myeloid_content_cells_sum<-data.frame(colSums( Xcell_Results_for_genes[which(rownames(Xcell_Results_for_genes) %in% myeloid_content_cells),] ))

lymphoid_myeloid_content_cells<-data.frame(lymphoid_content_cells_sum,myeloid_content_cells_sum)
colnames(lymphoid_myeloid_content_cells)<-c("lymphoid_content_cells_sum","myeloid_content_cells_sum")


write.table(lymphoid_myeloid_content_cells,"lymphoid_myeloid_content_cells",sep="\t")