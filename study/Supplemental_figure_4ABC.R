library(corrplot)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)  
library(data.table)
library(sleuth)
library(ComplexHeatmap)
library(splineTimeR)
#FUNCTION TO NEGATE %in%
"%nin%" = Negate("%in%")



#supplemental figure 4A
ggplot(Ordered_pheno, aes(x=Immune_subtypes, y=purity,fill = Immune_subtypes)) + 
  geom_violin(alpha = 0.5) +
  # ylim(0,14) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),axis.title.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_signif(comparisons = list(c("Immune cold", "Immune hot")), 
              map_signif_level=TRUE)


t.test(Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune hot",]$purity,Ordered_pheno[Ordered_pheno$Immune_subtypes %in% "Immune cold",]$purity)

# Supplemental figure 4B
set.seed(2021)
d <- Ordered_pheno
treatment <- d$Immune_subtypes
outcome <- d$purity  
original <- diff(tapply(outcome, treatment, mean))
mean(outcome[treatment=="Immune cold"])-mean(outcome[treatment=="Immune hot"])
#Permutation test
permutation.test <- function(treatment, outcome, n){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(original))/(n)
  return(list(result, distribution))
}
test1 <- permutation.test(treatment, outcome, 10000)
hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")
test1[[1]]


#Supplemental figure 4C lymphoid content  and myeloid content are calculated using Xcell R package and cell types were summed 
Corr_lymphoid_myeloid_purity<-data.frame(Ordered_pheno[,c("lymphoid content cells sum","myeloid content cells sum","purity")] )
colnames(Corr_lymphoid_myeloid_purity)<-c("lymphoid","myeloid","purity")
rownames(Corr_lymphoid_myeloid_purity)<-Ordered_pheno$CMO_ID
M <-cor(Corr_lymphoid_myeloid_purity)



corrplot(M, type="upper",method = 'circle', order = 'AOE',tl.col = "black",tl.cex = 0.8,
        col=colorRampPalette(c("blue","white","red"))(200))

