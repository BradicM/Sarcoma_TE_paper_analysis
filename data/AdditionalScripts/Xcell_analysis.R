library(data.table)
library(xCell)

#Perform  Xcell enrichment analysis 

exprMatrix = read.table("gene_expression_mtrx",header=TRUE,row.names=1, as.is=TRUE) #gene expression mtx output from RediscoverTE


exprMatrix<-exprMatrix[,-c(1)]

Xcell_Results<-xCellAnalysis(exprMatrix)

write.table(Xcell_Results,"Xcell_Results", sep="\t", quote=F) 
