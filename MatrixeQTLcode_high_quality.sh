###########################
#Final eQTL: 
#########################
R
setwd('/Volumes/Hyden/blh226/MatrixeQTL')
library(MatrixEQTL)
SNP<-read.table('SNP_filtered_157.txt', stringsAsFactors=FALSE, header=TRUE)



#SNP_filtered1<-names(SNP)

#SNP_filtered1<-SNP[rowSums(is.na(SNP))<80, ]
#write.table(SNP_filtered1, 'SNP_filtered_50%_HQ.txt', col.names=TRUE, row.names=TRUE)

#SNP_filtered2<-SNP[rowSums(is.na(SNP))<16, ]
#write.table(SNP_filtered2, 'SNP_filtered_90%_HQ.txt', col.names=TRUE, row.names=TRUE)



#row.names(SNP)<-SNP[,1]
#SNP<-SNP[,-1]


snpspos<-read.table('snpsloc_filtered_final.txt', stringsAsFactors=FALSE, header=TRUE)
#snpspos_subset<-snpspos[snpspos$snpid %in% rownames(SNP), ]
#write.table(snpspos_subset, 'snpsloc_filtered_final.txt', col.names=TRUE, row.names=FALSE) 



#set file locations and basic parameters
errorCovariance=numeric()
cisDist=1e6
base.dir='/Volumes/Hyden/blh226/MatrixeQTL'
useModel=modelANOVA
output_file_name_cis='Final Matrix eQTL HQ output'
output_file_name_tra='Final Matrix eQTL HQ output'


#Counts<-read.delim('/Volumes/Hyden/blh226/mapped_reads_final/FPKM/FPKM_combined_high_quality.txt', sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)
#FPKM<-read.delim('expression_final_filtered_5 HQ.txt', sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)
#filter<-match(names(Counts), names(FPKM))
#Expression<-FPKM[,filter]
#write.table(Expression, 'expression_final_filtered_5_HQ.txt', sep=" ", row.names=TRUE)
#expression<-read.table('expression_final_filtered_5_157.txt', stringsAsFactors=FALSE
#expression_20<-expression[apply(expression, 1, function(c)sum(c!=0))>20,]
#write.table(expression_20, 'expression_final_filtered_20_157.txt', row.names=TRUE, col.names=TRUE)

#genelocs<-read.table('genelocs_final_filtered_5.txt', stringsAsFactors=FALSE, header=TRUE)
#genelocs_subset_20<-genelocs[genelocs$Gene_ID %in% rownames(expression_20), ]
#write.table(genelocs_subset_20, 'genelocs_final_filtered_20.txt', col.names=TRUE, row.names=FALSE)


expression_file_name=paste(base.dir, '/expression_final_filtered_20_157.txt', sep="")
gene_location_file_name=paste(base.dir, '/genelocs_final_filtered_20.txt', sep="")



#snp<-read.table('SNP_filtered_final.txt', sep=" ", stringsAsFactors=FALSE, header=TRUE)
#names(snp)<-names(FPKM)[-181]
#filter<-match(names(Counts), names(snp))
#snpsnew<-snp[,filter]
#write.table(snpsnew, 'SNP_filtered_final_HQ.txt', sep=" ", row.names=TRUE)




SNP_file_name=paste(base.dir, "/SNP_filtered_157.txt", sep="")
snps_location_file_name=paste(base.dir, "/snpsloc_filtered_final.txt", sep="")
pvOutputThreshold_cis=2e-2
pvOutputThreshold_tra=1e-3

#cov<-read.table('covariate_final_withPCA.txt', sep=" ", stringsAsFactors=FALSE, header=TRUE)
#names(cov)<-names(FPKM)[-181]
#filter<-match(names(Counts), names(cov))
#covariates<-cov[,filter]
#write.table(covariates, 'covariate_final_withPCA_HQ.txt', sep=" ", row.names=FALSE)

covariates_file_name=paste(base.dir, '/covariate_final_withPCA_HQ_157.txt', sep="")

#read in files 
cvrt=SlicedData$new()
cvrt$fileDelimiter=' '
cvrt$fileOmitCharacters="NA"
cvrt$fileSkipRows=1
cvrt$fileSkipColumns=0
cvrt$LoadFile(covariates_file_name)

snps=SlicedData$new()
snps$fileDelimiter=' '
snps$fileOmitCharacters="NA"
snps$fileSkipRows=1
snps$fileSkipColumns=1
snps$fileSliceSize=2000
snps$LoadFile(SNP_file_name)

gene=SlicedData$new()
gene$fileDelimiter=' '
gene$fileOmitCharacters="NA"
gene$fileSkipRows=1
gene$fileSkipColumns=1
gene$fileSliceSize=2000
gene$LoadFile(expression_file_name)

genepos=read.table(gene_location_file_name, header=TRUE, stringsAsFactors=FALSE)
snpspos=read.table(snps_location_file_name, header=TRUE, stringsAsFactors=FALSE)

#filter for minor allele frequency 
maf.list = vector('list', length(snps))
for(i in 1:length(snps)) {
  slice = snps[[i]]
  maf.list[[i]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[i]] = pmin(maf.list[[i]],1-maf.list[[i]]);
}
maf = unlist(maf.list)

cat('SNPs before filtering:',nrow(snps))

#snps$RowReorderSimple(maf>0.05);
snps$RowReorder(maf>0.05);
cat('SNPs after filtering:',nrow(snps))

#normalize FPKM expression values 
for( i in 1:length(gene) ) {
  mat = gene[[i]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[i]] = mat;
}
rm(i, mat);

#run eQTL 
mefinal=Matrix_eQTL_main(snps=snps, gene=gene, cvrt=cvrt, output_file_name=output_file_name_tra, pvOutputThreshold=pvOutputThreshold_tra, useModel=useModel, errorCovariance=errorCovariance, verbose=TRUE, output_file_name.cis=output_file_name_cis, pvOutputThreshold.cis=pvOutputThreshold_cis, snpspos=snpspos, genepos=genepos, cisDist=cisDist, pvalue.hist='qqplot', min.pv.by.genesnp=TRUE, noFDRsaveMemory=FALSE)

#output results and figures 
write.table(mefinal$trans$eqtls, 'Final SNP coverage eQTL trans.txt')
write.table(mefinal$cis$eqtls, 'Final SNP coverage eQTL cis.txt')
write.table(mefinal$all$eqtls, 'Final SNP coverage all.txt')
pdf('mefinal_HQ.pdf', width=8, height=8)
plot(mefinal)
dev.off()
write.table(mefinal$trans$min.pv.snps, 'Final SNP coverage min_pv_snps_trans.txt')
write.table(mefinal$trans$min.pv.gene, 'Final SNP coverage min_pv_gene_trans.txt')
write.table(mefinal$cis$min.pv.snps, 'Final SNP coverage min_pv_snps_cis.txt')
write.table(mefinal$cis$min.pv.gene, 'Final SNP coverage min_pv_gene_cis.txt')

All_QTL<-rbind(mefinal$trans$eqtls, mefinal$cis$eqtls)
All_QTL_sig<-All_QTL[All_QTL$FDR<0.05,]
write.table(All_QTL_sig, 'Final eQTL SNP coverage HQ all.txt')

#All_QTL<-read.table('Final eQTL HQ all.txt', stringsAsFactors=FALSE)
#Wset<-grep('S20_', All_QTL$snps)
#Wset<-apply(All_QTL$snps, 1, grep, pattern='S20_')
#Wset<-grepl.sub(data=All_QTL, pattern='S20_', var=snps)
#Wsubset<-subset(All_QTL, select = grepl('S20_', All_QTL$snps))
Wsubset<-All_QTL_sig[grepl('S20_', All_QTL_sig$snps), ]
write.table(Wsubset, 'WeQTL.txt')

Zsubset<-All_QTL_sig[grepl('S15_', All_QTL_sig$snps), ]
write.table(Zsubset, 'ZeQTL.txt')

#subset for SDR
library(stringr)
two_columns<-str_split_fixed(Wsubset$snps, "_", n=2)
Wsubset2<-cbind(two_columns, Wsubset[,2:5])
WSDRsubset<-Wsubset2[as.numeric(as.character(Wsubset2[,2]))>2340000 & as.numeric(as.character(Wsubset2[,2]))<9070000,]
write.table(WSDRsubset, 'W-SDReQTL.txt')

two_columns<-str_split_fixed(Zsubset$snps, "_", n=2)
Zsubset2<-cbind(two_columns, Zsubset[,2:5])
ZSDRsubset<-Zsubset2[as.numeric(as.character(Zsubset2[,2]))>2340000 & as.numeric(as.character(Zsubset2[,2]))<6400000,]
write.table(ZSDRsubset, 'Z-SDReQTL.txt')

SDRsubset<-rbind(ZSDRsubset, WSDRsubset)
write.table(SDRsubset, 'SDReQTL.txt')

uniquegenes<-c()
holder<-SDRsubset$gene
number<-length(levels(as.factor(SDRsubset$gene)))
genes<-levels(as.factor(SDRsubset$gene))
for(i in 1:number){
  holder<-SDRsubset[SDRsubset$gene==genes[i],]
  uniquegenes<-rbind(uniquegenes, holder[which.min(holder$FDR),])
}
write.csv(uniquegenes, 'SDReQTL_uniquegenes.csv')

All_QTL_sig<-All_QTL[All_QTL$FDR<0.01,]
write.table(All_QTL_sig, 'Final QTL HQ FDR<0.01.txt')

All_QTL_sig<-All_QTL[All_QTL$FDR<0.05,]
write.table(All_QTL_sig_0.05, 'Final QTL HQ FDR<0.05.txt')

two_columns<-str_split_fixed(All_QTL_sig$snps, "_", n=2)
All_QTL_sig<-cbind(All_QTL_sig[,1], two_columns, All_QTL_sig[,2:5])
write.table(All_QTL_sig, 'Final eQTL SNP HQ all.txt')

trans_sig<-trans[trans$FDR<0.05,]
write.table(trans_sig, 'Final QTL trans FDR<0.05.txt')
two_columns<-str_split_fixed(trans_sig$snps, "_", n=2)
trans_sig<-cbind(trans_sig[,1], two_columns, trans_sig[,2:5])
write.table(trans_sig, 'Final eQTL SNP HQ trans.txt')



#Generate Plots
#loci<-levels(as.factor(All_QTL_sig[,1]))
loci<-levels(as.factor(trans[,1]))
#loci<-levels(as.factor(cis[,1]))

hits<-c()
for(i in 1:length(loci)){
	temp<-loci[i]
	value<-sum(trans[,1]==loci[i])
	hits<-c(hits, value)
	print(loci[i])
}


Summary<-cbind(loci, str_split_fixed(loci, "_", n=2), hits)
Summary<-as.data.frame(Summary)
names(Summary)<-c('Locus', 'Chr', 'Pos', 'eQTLs')
#need to sort summary by chr and locus
Chr<-str_split_fixed(Summary$Chr, "S", 2)
Summary$Chr<-Chr[,2]
Summary$Pos<-as.numeric(as.character(Summary$Pos))
Summary$Chr<-as.numeric(Summary$Chr)
Summary<-Summary[order(Summary$Pos),]
Summary<-Summary[order(Summary$Chr),]
write.table(Summary, 'eQTL_hits_per_locus_trans.txt', row.names=FALSE, col.names=TRUE)

#lineplot of summary
Summary$eQTLs<-as.numeric(as.character(Summary$eQTLs))


#lots of SNPs with extremely high numbers of QTLs
#Summary[Summary$eQTLs==max(Summary$eQTLs),]
#           Locus Chr   Pos eQTLs
#27383 S6_17605791   6 12043  3457
#27384 S6_17605813   6 12044  3457
#check SNP file
#SNP[rownames(SNP)=='S6_17605791',]
#SNP[rownames(SNP)=='S20_12462947',]
#these all have SNPs that only have a single example of a particular homozygote or heterozygote
#these are usually in the parents for some reason
#should exclude parents from eQTL mapping, would give much more consistant QTL results
#Summary[Summary$eQTLs>800,]
#87 loci
color<-rep(0, length(Summary$Chr))

color<-ifelse((Summary$Chr%%2==0)==TRUE, 'firebrick3', 'gold')
Summary$color<-as.character(color)
color<-ifelse((Summary$Chr==15 | Summary$Chr==20)==TRUE, 'deepskyblue', Summary$color)
Summary$color<-as.character(color)

pdf('Final eQTL Plot_trans.pdf', height=12, width=24)
barplot(Summary$eQTLs, xlab='Locus', ylab='No. of Associated Genes per locus', col=Summary$color, border=Summary$color)
axis(side=1, at=c(1000, 3500, 5800, 7630, 9550, 11800, 13550, 15500, 17200, 19050, 21500, 23050, 24500, 26150, 28200, 31100, 33900, 35300, 36500, 38000), labels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15Z', '16', '17', '18', '19', '15W'), tick=FALSE)
dev.off()

expression<-read.table('expression_final_filtered_20_157.txt', stringsAsFactors=FALSE
SNP[rownames(SNP)=='S234_41300',]
expression[rownames(expression)=='Sapur.019G046000.3.v5.1', ]

pdf('Top Cis-eQTL.pdf')
plot(jitter(as.numeric(SNP[rownames(SNP)=='S4_13836507',]), amount=.05), as.numeric(expression[rownames(expression)=='Sapur.004G144000.1.v5.1', ]), pch=19, main="S4_13836507 vs Sapur.004G144000.1", xlab="Genotype", ylab='Normalized Gene Expression', axes=FALSE)
box()
axis(side=2)
axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
reg<-lm(as.numeric(expression[rownames(expression)=='Sapur.004G144000.1.v5.1', ])~as.numeric(SNP[rownames(SNP)=='S4_13836507',]))
abline(reg, col='red')
dev.off()



top300<-read.table('Top_300_transeQTL.txt', header=TRUE)
#top300<-read.table('Top_300_transeQTL_diff_chr.txt', header=TRUE)
cov<-read.table('covariate_final_withPCA_HQ_157.txt', header=TRUE)
cov<-cov[5,]
for(i in 1:length(cov)){
	if(cov[i]==0){
		cov[i]='magenta'
	}else {
		cov[i]='blue'
	}
}

pdf('Top_64_eQTL_plots_diff_chr.pdf', width=25, height=13)
par(mar=c(2,4,1,1), mfrow=c(8,8))
for(i in 1:64){
	Name<-as.character(top300[i,1])
	Gene<-as.character(top300[i,4])
	plot(jitter(as.numeric(SNP[rownames(SNP)==Name,]), amount=.05), as.numeric(expression[rownames(expression)==Gene, ]), pch=19, main=paste(Name, 'vs', Gene), ylab='FPKM', axes=FALSE, col=adjustcolor(unlist(cov), 0.35))
	box()
	axis(side=2)
	axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
	reg<-lm(as.numeric(expression[rownames(expression)==Gene, ])~as.numeric(SNP[rownames(SNP)==Name,]))
	abline(reg, col='red')
}
dev.off()


##par(new=T)
####plot(as.factor(SNP[rownames(SNP)=='S4_13836507',]), as.numeric(expression[rownames(expression)=='Sapur.004G144000.1.v5.1', ]))

TF<-c()
holder<-All[All$gene==TFs$gene[2],]
holder[which.min(holder$FDR),]
for(i in 1:length(TFs$FDR)){
  holder<-All[All$gene==TFs$gene[i],]
  TF<-rbind(TF, holder[which.min(holder$FDR),])
}
write.csv(TF, 'TF.csv')


hotspots<-Summary[Summary$eQTLs>299,]
write.csv(hotspots, 'eQTL hotspots 300+ eQTLs.csv', col.names=TRUE, row.names=FALSE)

top_trans_qtl<-trans_sig[trans_sig[,1]=='S7_7498506',]
write.table(trans_sig[trans_sig[,1]=='S7_7498506',], 'S7_7498506 associated genes.txt', row.names=FALSE) 
pdf('Top_64_eQTL_plots_S7_7498506.pdf', width=25, height=13)
par(mar=c(2,4,1,1), mfrow=c(8,8))
for(i in 1:64){
	Name<-as.character(top_trans_qtl[i,1])
	Gene<-as.character(top_trans_qtl[i,4])
	plot(jitter(as.numeric(SNP[rownames(SNP)==Name,]), amount=.05), as.numeric(expression[rownames(expression)==Gene, ]), pch=19, main=paste(Name, 'vs', Gene), ylab='FPKM', axes=FALSE, col=adjustcolor(unlist(cov), 0.35))
	box()
	axis(side=2)
	axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
	reg<-lm(as.numeric(expression[rownames(expression)==Gene, ])~as.numeric(SNP[rownames(SNP)==Name,]))
	abline(reg, col='red')
}
dev.off()

second_top_trans_qtl<-trans_sig[trans_sig[,1]=='S19_1853770',]
write.table(second_top_trans_qtl, 'S19_1853770 associated genes.txt', row.names=FALSE) 
pdf('Top_64_eQTL_plots_S19_1853770.pdf', width=25, height=13)
par(mar=c(2,4,1,1), mfrow=c(8,8))
for(i in 1:64){
	Name<-as.character(second_top_trans_qtl[i,1])
	Gene<-as.character(second_top_trans_qtl[i,4])
	plot(jitter(as.numeric(SNP[rownames(SNP)==Name,]), amount=.05), as.numeric(expression[rownames(expression)==Gene, ]), pch=19, main=paste(Name, 'vs', Gene), ylab='FPKM', axes=FALSE, col=adjustcolor(unlist(cov), 0.35))
	box()
	axis(side=2)
	axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
	reg<-lm(as.numeric(expression[rownames(expression)==Gene, ])~as.numeric(SNP[rownames(SNP)==Name,]))
	abline(reg, col='red')
}
dev.off()

chr15Z<-All_QTL_sig[All_QTL_sig[,2]=='S15',]
chr15W<-All_QTL_sig[All_QTL_sig[,2]=='S20',]

chr15<-rbind(chr15Z, chr15W)
interesting_genes<-read.csv('/Volumes/Hyden/blh226/MatrixeQTL/SDR_eQTL_Interesting_genes.csv', stringsAsFactors=FALSE, header=TRUE)
set<-match(interesting_genes$ID, chr15$gene)
interesting_eQTL<-chr15[set,]

pdf('Interesting SDR eQTL.pdf', width=50, height=26)
par(mar=c(2,4,1,1), mfrow=c(12,12))
for(i in 1:138){
	Name<-as.character(interesting_eQTL[i,1])
	Gene<-as.character(interesting_eQTL[i,4])
	plot(jitter(as.numeric(SNP[rownames(SNP)==Name,]), amount=.05), as.numeric(expression[rownames(expression)==Gene, ]), pch=19, main=paste(Name, 'vs', Gene), ylab='FPKM', axes=FALSE, col=adjustcolor(unlist(cov), 0.35))
	box()
	axis(side=2)
	axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
	reg<-lm(as.numeric(expression[rownames(expression)==Gene, ])~as.numeric(SNP[rownames(SNP)==Name,]))
	abline(reg, col='red')
}

dev.off()

X4G1044<-All_QTL_sig[All_QTL_sig[,4]=='Sapur.004G104400.1.v5.1',]
write.table(X4G1044, 'Sapur.004G104400.1 eQTL.txt')
#edit in text editor
GENE<-read.table('Sapur.004G104400.1 eQTL.txt', stringsAsFactors=FALSE)
pdf('Sapur.004G104400.1 eQTL.pdf', width=12, height=6)
par(mar=c(2,4,1,1), mfrow=c(3,3))
for(i in 1:5){
	Name<-as.character(GENE[i,1])
	Gene<-as.character(GENE[i,4])
	plot(jitter(as.numeric(SNP[rownames(SNP)==Name,]), amount=.05), as.numeric(expression[rownames(expression)==Gene, ]), pch=19, main=paste(Name, 'vs', Gene), ylab='FPKM', axes=FALSE, col=adjustcolor(unlist(cov), 0.35))
	box()
	axis(side=2)
	axis(side=1, at=c(0,1,2), labels=c('AA', 'AB', 'BB'))
	reg<-lm(as.numeric(expression[rownames(expression)==Gene, ])~as.numeric(SNP[rownames(SNP)==Name,]))
	abline(reg, col='red')
}
dev.off()

#calculating dominance deviations
GENE<-"Sapur.001G000400.3.v5.1"
LOCUS<-"S15_5678743"
table<-cbind(t(SNP[rownames(SNP)==LOCUS,]), t(expression[rownames(expression)==GENE,]))
table<-as.data.frame(table)
#table[,1]<-as.factor(table[,1])
reg<-lm(as.numeric(table[,2])~as.numeric(table[,1]))
mu<-as.numeric(coef(reg)[1] + coef(reg)[2])
Gab<-mean(table[table[,1]==1,2], na.rm=TRUE)
Gaa<-mean(table[table[,1]==0,2], na.rm=TRUE)
Gbb<-mean(table[table[,1]==2,2], na.rm=TRUE)
a1<-as.numeric(coef(reg)[1])
a2<-as.numeric(coef(reg)[1] + 2*(coef(reg)[2]))
#delta11<-Gaa-a1
delta12<-Gab-mu
#delta22<-Gbb-a2
pval<-pnorm(mu, mean=Gab, sd(table[table[,1]==1,2], na.rm=TRUE))


interesting_eQTL$het_dom_dev<-rep(NA, length(interesting_genes[,1]))
interesting_eQTL$dom_dev_percent<-rep(NA, length(interesting_genes[,1]))
interesting_eQTL$pval<-rep(NA, length(interesting_genes[,1]))
for(i in 1:length(interesting_eQTL[,1])){
	GENE<-interesting_eQTL[i,4]
	LOCUS<-interesting_eQTL[i,1]
	table<-cbind(t(SNP[rownames(SNP)==LOCUS,]), t(expression[rownames(expression)==GENE,]))
	table<-as.data.frame(table)
	reg<-lm(as.numeric(table[,2])~as.numeric(table[,1]))
	mu<-as.numeric(coef(reg)[1] + coef(reg)[2])
	Gab<-mean(table[table[,1]==1,2], na.rm=TRUE)
	delta12<-Gab-mu
	percent<-delta12/mu
	pval<-pnorm(mu, mean=Gab, sd(table[table[,1]==1,2], na.rm=TRUE))
	if(pval > 0.5){
		pval<-1-pval
	}
	interesting_eQTL$het_dom_dev[i]<-delta12
	interesting_eQTL$dom_dev_percent[i]<-percent
	interesting_eQTL$pval[i]<-pval
}
write.csv(interesting_eQTL, 'Interesting genes with dominance deviations.csv', row.names=FALSE)


#Circos plot
names(chr15)<-c('SNP', 'chr', 'pos', 'gene', 'statistic', 'pvalue', 'FDR')
Chr<-str_split_fixed(chr15$chr, "S", 2)
chr15$chr<-as.numeric(Chr[,2])
GFF3<-read.table('Spurpurea_519_v5.1.gene.gff3', header=FALSE, stringsAsFactors=FALSE)
GFF3<-GFF3[GFF3$V3=='mRNA', ]
chromosomes<-str_split_fixed(GFF3$V1, "Chr", 2)
chromosomes<-chromosomes[,2]

genes<-str_split_fixed(GFF3$V9, "=", 2)
genes<-str_split_fixed(genes[,2], ";", 2)
genes<-genes[,1]
geneinfo<-cbind(chromosomes, GFF3$V4, genes)
geneinfo[,2]<-as.numeric(geneinfo[,2])
write.table(geneinfo, 'geneinfo.txt', row.names=FALSE, col.names=TRUE)
matchy<-match(chr15[,4], geneinfo[,3])
chr15info<-cbind(chr15, geneinfo[matchy,])
for(i in 1:length(chr15info[,2])){
	if(chr15info[i,2]==15){
	chr15info[i,2]='15Z'
	}
}
for(i in 1:length(chr15info[,2])){
	if(chr15info[i,2]==20){
	chr15info[i,2]='15W'
	}
}
chr15info$chromosomes<-as.character(chr15info$chromosomes)
for(i in 1:length(chr15info[,8])){
	if(is.na(chr15info[i,8])){
		chr15info[i,8]="T"
	}
}

chr_table<-read.csv('chr_lengths.csv', stringsAsFactors=FALSE, header=TRUE)
chr15info$pos<-as.numeric(as.character(chr15info$pos))
chr15info$V2<-as.numeric(as.character(chr15info$V2))
SDRinfo<-chr15info[chr15info$pos>2340000,]
Zinfo<-SDRinfo[SDRinfo$chr=='15Z',]
Zinfo<-Zinfo[Zinfo$pos<6700000,]
Winfo<-SDRinfo[SDRinfo$chr=='15W',]
Winfo<-Winfo[Winfo$pos<9070000,]
SDRinfo<-rbind(Winfo, Zinfo)


nuc1<-cbind(SDRinfo$chr, as.numeric(as.character(SDRinfo$pos)), as.numeric(as.character(SDRinfo$pos))+1)
nuc2<-cbind(SDRinfo$chromosomes, as.numeric(as.character(SDRinfo$V2)), as.numeric(as.character(SDRinfo$V2))+1)
nuc1<-as.data.frame(nuc1)
nuc2<-as.data.frame(nuc2)
nuc1[,4]<-rep(0, length(nuc1[,1]))
for(i in 1:length(nuc1[,1])){
	if(nuc1[i,1]=='15Z'){
		nuc1[i,4]='blue'
	}else{
		nuc1[i,4]='magenta'
	}
}


library(circlize)
library(migest)
library(dplyr)
library(scales)

n=nrow(SDRinfo)
par(mar=rep(0,4))
circos.clear()
col_text <- "black"
circos.par("track.height"=0.8,gap.degree=1,cell.padding=c(0,0,0,0))
## Basic circos graphic parameters
#circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)
circos.initialize(factors = chr_table$Chr, xlim = cbind(chr_table$Start, chr_table$End))
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim),mean(ylim),chr,cex=1,col=col_text,
facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)

brk <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)*10^6
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)


for(i in 1:length(nuc1[,1])){
	circos.link(sector.index1=as.character(nuc1[i,1]), point1=as.numeric(as.character(nuc1[i,2])), sector.index2=as.character(nuc2[i,1]), point2=as.numeric(as.character(nuc2[i,2])), col=alpha(nuc1[i,4], 0.35), lwd=0.1)
}







