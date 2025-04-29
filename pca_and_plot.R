# R script used to perform and plot a PCA from the covariance matrix obtained with PCAngsd

# Author: Carlos Leiva

#Load the libraries
library(RcppCNPy)
library(tibble)
library(ggplot2)


#Load the covariance matrix
cov<-as.matrix(read.table("Aca_198_ind_thin10kSNPs_PCAngsd.cov"))


pop <- c(rep("Australia", 19), rep("French Polynesia", 13), rep("Fiji", 4), rep("Hawaii", 9), rep("Japan", 26),
         rep("Marshall", 7), rep("French Polynesia", 13), rep("Japan", 41), rep("PNG", 4), rep("Philippines", 4), rep("Pohnpei", 9),
         rep("French Polynesia", 23), rep("Vanuatu", 16), rep("Vietnam", 10))


mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4


ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) +
  geom_point(size=4) +
  theme_classic() +
  theme(#legend.position = c(0.8, 0.8), 
    legend.background = element_rect(color = 'black'),
    legend.title = element_blank()) +
  xlab(paste0('PC1 (', sprintf('%0.1f%% explained var.', 100*mme.pca$values[1]/pca.eigenval.sum), ")"))+
  ylab(paste0('PC2 (', sprintf('%0.1f%% explained var.', 100*mme.pca$values[2]/pca.eigenval.sum), ")"))

