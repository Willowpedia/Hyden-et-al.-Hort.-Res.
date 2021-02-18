##########################################
#Calculate FPKM 
###########################################
#get rid of columns and rows with names 
counts<-allfilescombined[,-1]
counts1<-mapply(counts, FUN=as.numeric)
counts2<-matrix(data=counts1, ncol=183, nrow=55612)
genes<-as.integer(as.character(allfiles[[1]]$Length))
genes1<-cbind(allfiles[[1]]$Geneid, genes)
colnames(genes1)<-c('Gene_ID', 'Length')

#unlist counts and turn into matrix 
library(edgeR)
myDGEList<-DGEList(counts=counts2, genes=genes1)
newDGEList<-calcNormFactors(myDGEList)
fpkmMatrix<-rpkm(newDGEList)
setwd('/Volumes/Hyden/blh226/mapped_reads_final/FPKM')
rownames(fpkmMatrix)<-allfiles[[1]]$Geneid
write.table(fpkmMatrix, 'FPKM_combined.txt', col.names=TRUE, row.names=TRUE)
write.table(fpkmMatrix, 'FPKM_combined_nonames.txt', col.names=FALSE, row.names=FALSE)

