# R script used to calculate and plot Pairwise Fst distances among populations.

# Author: Carlos Leiva

################
# PAIRWISE FST #
################

library(adegenet)
library(vcfR)
library(tidyverse)
library(parallel)
library(ggplot2)
library(StAMPP)
library(pheatmap)

#Read the vcf and transform to genlight

snps <- read.vcfR("Aca_198_ind_think10kSNPs.recode.vcf")
snps_gl<- vcfR2genlight(snps)

#Add pop info

pop <- c(rep("Australia", 19), rep("BoraBora", 13), rep("Fiji", 4), rep("Hawaii", 9), rep("Kagoshima", 26),
         rep("Marshall", 7), rep("Moorea", 13), rep("Okinawa", 41), rep("PNG", 4), rep("Philippines", 4), rep("Pohnpei", 9),
         rep("Raiatea",8), rep("Tahiti", 15), rep("Vanuatu", 16), rep("Vietnam", 10))

snps_gl@pop <- as.factor(pop)


#run the pairwise fst analysis with 1000 bootstrap replicates.

fsttable_pval<-stamppFst(snps_gl, nclusters = 4)#This is with pvalues and bootstraps. It did pval, 100 is default!

fsttable_pval<-stamppFst(snps_gl, nboots = 1000, percent = 95, nclusters = 4)#This is with pvalues and bootstraps.

fsttable_pval

write.csv(fsttable_pval$Fsts, "pairwise_fst.csv")
write.csv(fsttable_pval$Pvalues, "pairwise_fst_pval.csv")
write.csv(fsttable_pval$Bootstraps, "pairwise_fst_bootstraps.csv")


#Plot fsts

FSTtab<-fsttable_pval$Fsts

diag(FSTtab)<-0

f <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

FSTtabcompl<-f(M)

M <- FSTtab
for(i in 1:nrow(M)) {for(j in 1:i) {M[i,j]=M[j,i] }}
M

m[upper.tri(m)] <- t(m)[upper.tri(m)]

pheatmap(FSTtabcompl)
