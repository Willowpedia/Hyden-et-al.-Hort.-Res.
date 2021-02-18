#test alignment 
~/programs/ShortStack-3.8.3/ShortStack --bowtie_cores 24 --sort_mem 54000M --dicermax 25 --mincov 12 --readfile fastq/11993.1.231799.ACAGT.filter-SMRNA.fastq.gz --genomefile /Volumes/Hyden/blh226/Spurpurea/v5.1/assembly/Spurpurea_519_v5.0.fa 
#align and analyze all files 
~/programs/ShortStack-3.8.3/ShortStack --bowtie_cores 24 --sort_mem 54000M --dicermax 25 --mincov 12 --readfile fastq/*.filter-SMRNA.fastq.gz --genomefile /Volumes/Hyden/blh226/Spurpurea/v5.1/assembly/Spurpurea_519_v5.0.fa 
#later work 
#create fasta of small rnas 
#look at rfam annotations, if number of miRNA loci is small, can try manually searching and annotating as well.  in R extract column of major sequence and search databases.  
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm
#do an eQTL analysis? 
#use data in network analysis? 
#analysis of data
#RPKM
R
#read in locus names
ShortStack_All<-read.table("/Volumes/Hyden/blh226/smrna/f2_smrna_catkin_analysis/ShortStack_1562789001/ShortStack_All.gff3")
#read in counts
counts<-read.table("/Volumes/Hyden/blh226/smrna/f2_smrna_catkin_analysis/ShortStack_1562789001/Counts.txt", stringsAsFactors=FALSE)
#get lengths of loci 
ShortStack_All[,10]<-abs(ShortStack_All[,5]-ShortStack_All[,4])
#get rid of unmapped counts, columns and rows with names 
counts<-counts[-1,]
counts<-counts[,-1]
counts<-counts[,-1]
counts<-counts[1:266272,]
#unlist counts and turn into matrix 
mat<-matrix(NA,nrow=266272,ncol=45)
for(i in 1:ncol(counts)){
	mat[,i]<-as.numeric(unlist(counts[,i]))
	}
#calculate RPKM values
library(edgeR)
myDGEList<-DGEList(counts=mat, genes=ShortStack_All)
newDGEList<-calcNormFactors(myDGEList)
rpkmMatrix<-rpkm(newDGEList)
write.table(rpkmMatrix, 'RPKM.txt')
metadata<-read.csv('/Volumes/Hyden/blh226/smrna/smRNA-metadata-2019-07-08-CHC.csv' , stringsAsFactors=FALSE)
#rename genotypes to be consistant with other files 
Genotype<-metadata[2:46,10]
for(i in 1:length(Genotype)){
 Genotype[i]<-paste0('X',Genotype[i])
}
RPKMnew<-rbind(Genotype, RPKM)
write.table(RPKMnew, 'RPKM.new.txt')
#create sex covariate file 
Sex<-c(0)
for(i in 1:length(Gender)){
 if(Gender[i]=='Female'){
 Sex<-cbind(Sex, 0)
 }
 else if(Gender[i]=='Male'){
 Sex<-cbind(Sex, 1)
 }
 else{
 cbind(Sex, NA)
 }
}
Covariate<-rbind(Genotype, Sex)
write.table(Covariate, 'Covariate.txt')
#remove genotypes from RPKM file that do not have SNP files that are ready (94001, 94003, 94006)
RPKMnew<-RPKMnew[,-16]
RPKMnew<-RPKMnew[,-14]
RPKMnew<-RPKMnew[,-5]
write.table(RPKMnew, 'RPKM.new.subset.txt')
#prepare SNP file 
#read in file with all samples 
hmp<-read.table('/Volumes/HAL9000/genomics/salix_ngs_fastq/fastq_gbs/Hyden/MatrixeQTL_files/15ZW_eco_ape_ABH_recode.hmp.txt', header=TRUE, stringsAsFactors=FALSE)
#remove info columns, recode AHB to 012
SNP<-cbind(hmp[,1], hmp[,12:ncol(hmp)])
SNP[SNP=="A"]<-0
SNP[SNP=="H"]<-1
SNP[SNP=="B"]<-2
write.table(SNP, file="SNP.txt", col.names=TRUE, row.names=FALSE)
SNPsubset2<-as.character(SNP[,1])
for(i in 1:length(Genotypes)){
	SNPsubset2<-cbind(SNPsubset2, SNP[,names(SNP)%in%Genotypes[i]])
}
write.table(SNPsubset2, 'SNP.subset.txt', row.names=FALSE, col.names=FALSE)
#create genelocs table
genelocs<-read.table('/Volumes/Hyden/blh226/smrna/f2_smrna_catkin_analysis/ShortStack_1562789001/ShortStack_All.gff3', stringsAsFactors=FALSE)
genelocs<-cbind(genelocs[,9], genelocs[,1], genelocs[,4], genelocs[,5])
names(genelocs)<-c('Gene', 'chr', 's1', 's2')
write.table(genelocs, 'geneloc.txt', row.names=FALSE, col.names=TRUE)

#Matrix eQTL 
library(MatrixEQTL)
useModel=modelLINEAR_CROSS
SNP_file_name='SNP.subset.txt'
expression_file_name='RPKM.new.subset.txt'
covariates_file_name='Covariate.txt'
snps_location_file_name='snpsloc.txt'
gene_location_file_name='geneloc.txt'
output_file_name_cis='cis_output'
output_file_name_tra='trans_output'
pvOutputThreshold_tra=0.01
pvOutputThreshold_cis=0.02
errorCovariance=numeric()
cisDist=1e6

snps=SlicedData$new()
snps$fileDelimiter=' '
snps$fileOmitCharacters='NA'
snps$fileSkipRows=1
snps$fileSkipColumns=1
snps$fileSliceSize=2000
snps$LoadFile(SNP_file_name)

gene=SlicedData$new()
gene$fileDelimiter=' '
gene$fileOmitCharacters='NA'
gene$fileSkipRows=1
gene$fileSkipColumns=1
gene$fileSliceSize=2000
gene$LoadFile(expression_file_name)

cvrt=SlicedData$new()
cvrt$fileDelimiter=' '
cvrt$fileOmitCharacters='NA'
cvrt$fileSkipRows=1
cvrt$fileSkipColumns=0
cvrt$LoadFile(covariates_file_name)

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)


me=Matrix_eQTL_main(snps = snps, gene=gene, cvrt=cvrt, output_file_name=output_file_name_tra, pvOutputThreshold=pvOutputThreshold_tra, useModel=useModel, errorCovariance=errorCovariance, verbose=TRUE, output_file_name.cis=output_file_name_cis, pvOutputThreshold.cis=pvOutputThreshold_cis, snpspos=snpspos, genepos=genepos, cisDist=cisDist, pvalue.hist='qqplot', min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE)
#Ran out of memory, try setting no FDR to true to save memory, also try ANOVA as it is not entirely clear what Linear cross measures 
useModel=modelANOVA
me<-Matrix_eQTL_main(snps = snps, gene=gene, cvrt=cvrt, output_file_name=output_file_name_tra, pvOutputThreshold=pvOutputThreshold_tra, useModel=useModel, errorCovariance=errorCovariance, verbose=TRUE, output_file_name.cis=output_file_name_cis, pvOutputThreshold.cis=pvOutputThreshold_cis, snpspos=snpspos, genepos=genepos, cisDist=cisDist, pvalue.hist='qqplot', min.pv.by.genesnp=FALSE, noFDRsaveMemory=TRUE)

unlink(output_file_name_tra)
unlink(output_file_name_cis)
me$time.in.sec
me$cis
me$trans
plot(me)
#no cis effects, many highly significant trans effects, indicates that more covariates are necessary 
me$all

#PCA analysis on RPKM
RPKM<-read.table('RPKM.new.subset.txt', stringsAsFactors=FALSE)
t_RPKM<-t(RPKM)
t_RPKM<-t_RPKM[-1,-1]
t_RPKM<-apply(t_RPKM, 2, as.numeric)
t_RPKM<-t_RPKM[,apply(t_RPKM, 2, var, na.rm=TRUE) !=0]
write.table(t_RPKM, 't_RPKM.txt')
myPCA<-prcomp(t_RPKM, center=TRUE, scale.=TRUE)
summary(myPCA)
#first PCA explains 13% of variance
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(myPCA)
#install.packages('factoextra')
library(factoextra)
fviz_eig(myPCA)
write.table(t(myPCA$x[,1]), 'PCA1.txt', col.names=FALSE, row.names=FALSE)
write.table(t(myPCA$x[,2]), 'PCA2.txt', col.names=FALSE, row.names=FALSE)
write.table(t(myPCA$x[,3]), 'PCA3.txt', col.names=FALSE, row.names=FALSE)
write.table(t(myPCA$x[,4]), 'PCA4.txt', col.names=FALSE, row.names=FALSE)
#paste PCA 1-4 in covariates file and re-do matrix eqtl 
#still too many significant p-values, try filtering out regions with low coverage or few samples with genes mapping, maybe only try miRNA subset 
#also cannot find locations of eqtls, output file names for cis and trans should have the same names 

