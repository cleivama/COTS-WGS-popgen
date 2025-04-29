# R script used to plot dendrogram and heatmap from relatedness data from vcftools --relatedness2

# Author: Carlos Leiva

library(data.table)
library(reshape2)
library(ggdendro)
library(pheatmap)


relations <- data.frame(fread('out.relatedness2'))

relate <- relations[,c("INDV1","INDV2","RELATEDNESS_PHI")]
rel <- dcast(relate, INDV1~INDV2)
row.names(rel)<-rel$INDV1
rel <- rel[,-1]

d <- dist(as.matrix(rel))
hc <- hclust(d)
ggdendrogram(hc, size=12)

pheatmap(rel)
