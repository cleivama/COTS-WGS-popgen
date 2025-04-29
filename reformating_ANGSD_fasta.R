# R script used to transform haplotype calls from ANGSD (-doHaploCall 2) into a fasta file to perform phylogenetic analyses with iqtree2.

# Authors: Carlos Leiva, based on a previous script from Héctor Torrado

#Read haplo file and remove the first three columns, which we don't want on our alignment
dt<-read.table("Aca_ALLhaplo_allSNPs.haplo", header = T)[, -c(1,2,3)]


#Calculate the missing data percentage per row (Position)
library(dplyr)
dt<-dt %>%
  mutate(missing_perc = rowMeans(select(., 1: length(colnames(dt)))=="N") * 100)


#Filter for 50% presence
dtf<- dt[dt$missing_perc<50,]
#Filter for 30% presence
dtf30<- dt[dt$missing_perc<70,]
#If the file is big, is good for your RAM to remove from memory the original file and clean a bit after the first filter
 rm(dt)
 gc() 


#Filters for 80%, 90% and 100% presence
dtf80<- dtf[dtf$missing_perc<20, -length(colnames(dtf))]
dtf90<- dtf[dtf$missing_perc<10, -length(colnames(dtf))]
dtf100<- dtf[dtf$missing_perc==0, -length(colnames(dtf))]


#Read a file with the names for the individuals. Can be just the bamlist or something more clean and easy to read (preferred)
nomz<-read.table("nomz_addCA_NoCont_addINDIAN.txt")


#Transpose the datasets
Tdtf<-t(dtf[ ,-length(colnames(dtf))])
Tdtf30<-t(dtf30[ ,-length(colnames(dtf))])
Tdtf80<-t(dtf80)
Tdtf90<-t(dtf90)
Tdtf100<-t(dtf100)


#Add names and format to look like a fasta alignment
Tdtf302<-cbind(paste(">",nomz$V1,"\n", sep=""), Tdtf30)
Tdtf2<-cbind(paste(">",nomz$V1,"\n", sep=""), Tdtf)
Tdtf802<-cbind(paste(">",nomz$V1,"\n", sep=""), Tdtf80)
Tdtf902<-cbind(paste(">",nomz$V1,"\n", sep=""), Tdtf90)
Tdtf1002<-cbind(paste(">",nomz$V1,"\n", sep=""), Tdtf100)


#Export as a file
write.table(Tdtf302, "T_mod_ANGSD_Haplo_03filt.fasta", col.names = F, row.names = F, sep = "", quote = F )
write.table(Tdtf2, "T_mod_ANGSD_Haplo_05filt.fasta", col.names = F, row.names = F, sep = "", quote = F )
write.table(Tdtf802, "T_mod_ANGSD_Haplo_08filt.fasta", col.names = F, row.names = F, sep = "", quote = F )
write.table(Tdtf902, "T_mod_ANGSD_Haplo_09filt.fasta", col.names = F, row.names = F, sep = "", quote = F )
write.table(Tdtf1002, "T_mod_ANGSD_Haplo_1filt.fasta", col.names = F, row.names = F, sep = "", quote = F )


