cd /Volumes/HAL9001/Hyden/bisulfite-seq/SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results
pigz -d *.cov.gz
cd ../../../
#test
#DMRfinder-master/combine_CPG_sites.py -o test.csv SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/HASZB.read1_bismark_bt2_pe.deduplicated.bismark.cov SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/HASZC.read1_bismark_bt2_pe.deduplicated.bismark.cov


#DMRfinder-master/combine_CPG_sites.py -o results.txt SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/*.cov
#doesn't work, uses too much memory
#try low memory option
for file in *.bismark.cov
do
	baseFilename=`basename $file .bismark.cov`
	echo "${baseFilename}"
	sort -k 1n "${baseFilename}".bismark.cov > "${baseFilename}".bismark.sorted.cov
done

#test
#DMRfinder-master/combine_CPG_sites.py -v -o testresults.txt SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/HASZB.read1_bismark_bt2_pe.deduplicated.bismark.sorted.cov SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/HASZC.read1_bismark_bt2_pe.deduplicated.bismark.sorted.cov -b -e chr_order.txt 


DMRfinder-master/combine_CPG_sites.py -v -o results.txt SalpurSequencing_6_FD/Bisulfite-seq_Analysis/Results/*.cov -b -e chr_order.txt

Rscript DMRfinder-master/findDMRs.r -v -i results.txt -o DMR.txt -n Male,Female HATCU,HATCW,HATCX,HATCY,HATCZ,HATGA,HATGB,HATGC,HATGG,HATGH,HATGN,HATGO,HATGP,HATGS,HATGT,HATGU HASZB,HASZC,HASZG,HASZH,HASZN,HASZO,HASZP,HASZS,HASZT,HASZU,HASZW,HASZX,HASZY,HASZZ,HATAA,HATAB

Rscript DMRfinder-master/findDMRs.r -v -i results.txt -o DMR_diversity.txt -n Male,Female HATCU,HATCW,HATCX,HATCY,HATCZ,HATGA,HATGB,HATGC,HATGG,HATGH HASZB,HASZC,HASZG,HASZH,HASZN,HASZO,HASZP,HASZS,HASZT,HASZU

Rscript DMRfinder-master/findDMRs.r -v -i results.txt -o DMR_F2.txt -n Male,Female HATGN,HATGO,HATGP,HATGS,HATGT,HATGU HASZW,HASZX,HASZY,HASZZ,HATAA,HATAB

Rscript DMRfinder-master/findDMRs.r -v -i results.txt -o DMR_H_vs_F.txt -n Monoecious,Female HATGW HASZB,HASZC,HASZG,HASZH,HASZN,HASZO,HASZP,HASZS,HASZT,HASZU,HASZW,HASZX,HASZY,HASZZ,HATAA,HATAB

Rscript DMRfinder-master/findDMRs.r -v -i results.txt -o DMR_H_vs_M.txt -n Male,Female HATCU,HATCW,HATCX,HATCY,HATCZ,HATGA,HATGB,HATGC,HATGG,HATGH,HATGN,HATGO,HATGP,HATGS,HATGT,HATGU HATGW

#match loci with nearest genes
R
library(MALDIquant)
gff3<-read.table('Spurpurea_519_v5.1.gene.gff3', stringsAsFactors=FALSE, header=FALSE)
gff3<-gff3[gff3[,3]=='gene',]
DMR<-read.table('DMR.txt', stringsAsFactors=FALSE, header=TRUE)
DMRgenes<-c()
for(i in 1:length(DMR$chr)){
	print(i)
	Chr<-gff3[gff3[,1]==DMR[i,1],]
	x<-DMR[i,3]
	y<-DMR[i,2]
	pos<-Chr[Chr[,7]=='+',]
	pos<-pos[pos[,4] > x,]
	pos<-pos[order(pos[,4]),]
	neg<-Chr[Chr[,7]=='-',]
	neg<-neg[neg[,5] <y, ]
	neg<-neg[order(neg[,5]),]
	sense<-pos[match.closest(x, pos[,4]),]
	antisense<-neg[match.closest(y, neg[,5]),]
	if(is.na(sense[1,4])){
		both<-cbind(DMR[i,], antisense)
	}else if(sense[,4]-x < y-antisense[,5] | is.na(antisense[1,5])){
		both<-cbind(DMR[i,], sense)
	}else{
		both<-cbind(DMR[i,], antisense)
	}
	DMRgenes<-rbind(DMRgenes, both)
}
write.csv(DMRgenes, 'DMR_nearest_genes.csv', row.names=FALSE)
write.table(DMRgenes, 'DMR_nearest_genes.txt', row.names=FALSE)

DMR<-read.table('DMR_diversity.txt', stringsAsFactors=FALSE, header=TRUE)
DMRgenes<-c()
for(i in 1:length(DMR$chr)){
	print(i)
	Chr<-gff3[gff3[,1]==DMR[i,1],]
	x<-DMR[i,3]
	y<-DMR[i,2]
	pos<-Chr[Chr[,7]=='+',]
	pos<-pos[pos[,4] > x,]
	pos<-pos[order(pos[,4]),]
	neg<-Chr[Chr[,7]=='-',]
	neg<-neg[neg[,5] <y, ]
	neg<-neg[order(neg[,5]),]
	sense<-pos[match.closest(x, pos[,4]),]
	antisense<-neg[match.closest(y, neg[,5]),]
	if(is.na(sense[1,4])){
		both<-cbind(DMR[i,], antisense)
	}else if(sense[,4]-x < y-antisense[,5] | is.na(antisense[1,5])){
		both<-cbind(DMR[i,], sense)
	}else{
		both<-cbind(DMR[i,], antisense)
	}
	DMRgenes<-rbind(DMRgenes, both)
}
write.csv(DMRgenes, 'DMR_diversity_nearest_genes.csv', row.names=FALSE)
write.table(DMRgenes, 'DMR_diversity_nearest_genes.txt', row.names=FALSE)

DMR<-read.table('DMR_F2.txt', stringsAsFactors=FALSE, header=TRUE)
DMRgenes<-c()
for(i in 1:length(DMR$chr)){
	print(i)
	Chr<-gff3[gff3[,1]==DMR[i,1],]
	x<-DMR[i,3]
	y<-DMR[i,2]
	pos<-Chr[Chr[,7]=='+',]
	pos<-pos[pos[,4] > x,]
	pos<-pos[order(pos[,4]),]
	neg<-Chr[Chr[,7]=='-',]
	neg<-neg[neg[,5] <y, ]
	neg<-neg[order(neg[,5]),]
	sense<-pos[match.closest(x, pos[,4]),]
	antisense<-neg[match.closest(y, neg[,5]),]
	if(is.na(sense[1,4])){
		both<-cbind(DMR[i,], antisense)
	}else if(sense[,4]-x < y-antisense[,5] | is.na(antisense[1,5])){
		both<-cbind(DMR[i,], sense)
	}else{
		both<-cbind(DMR[i,], antisense)
	}
	DMRgenes<-rbind(DMRgenes, both)
}
write.csv(DMRgenes, 'DMR_F2_nearest_genes.csv', row.names=FALSE)
write.table(DMRgenes, 'DMR_F2_nearest_genes.txt', row.names=FALSE)

DMR<-read.table('DMR_H_vs_F.txt', stringsAsFactors=FALSE, header=TRUE)
DMRgenes<-c()
for(i in 1:length(DMR$chr)){
	print(i)
	Chr<-gff3[gff3[,1]==DMR[i,1],]
	x<-DMR[i,3]
	y<-DMR[i,2]
	pos<-Chr[Chr[,7]=='+',]
	pos<-pos[pos[,4] > x,]
	pos<-pos[order(pos[,4]),]
	neg<-Chr[Chr[,7]=='-',]
	neg<-neg[neg[,5] <y, ]
	neg<-neg[order(neg[,5]),]
	sense<-pos[match.closest(x, pos[,4]),]
	antisense<-neg[match.closest(y, neg[,5]),]
	if(is.na(sense[1,4])){
		both<-cbind(DMR[i,], antisense)
	}else if(sense[,4]-x < y-antisense[,5] | is.na(antisense[1,5])){
		both<-cbind(DMR[i,], sense)
	}else{
		both<-cbind(DMR[i,], antisense)
	}
	DMRgenes<-rbind(DMRgenes, both)
}
write.csv(DMRgenes, 'DMR_H_vs_F_nearest_genes.csv', row.names=FALSE)
write.table(DMRgenes, 'DMR_H_vs_F_nearest_genes.txt', row.names=FALSE)

DMR<-read.table('DMR_H_vs_M.txt', stringsAsFactors=FALSE, header=TRUE)
DMRgenes<-c()
for(i in 1:length(DMR$chr)){
	print(i)
	Chr<-gff3[gff3[,1]==DMR[i,1],]
	x<-DMR[i,3]
	y<-DMR[i,2]
	pos<-Chr[Chr[,7]=='+',]
	pos<-pos[pos[,4] > x,]
	pos<-pos[order(pos[,4]),]
	neg<-Chr[Chr[,7]=='-',]
	neg<-neg[neg[,5] <y, ]
	neg<-neg[order(neg[,5]),]
	sense<-pos[match.closest(x, pos[,4]),]
	antisense<-neg[match.closest(y, neg[,5]),]
	if(is.na(sense[1,4])){
		both<-cbind(DMR[i,], antisense)
	}else if(sense[,4]-x < y-antisense[,5] | is.na(antisense[1,5])){
		both<-cbind(DMR[i,], sense)
	}else{
		both<-cbind(DMR[i,], antisense)
	}
	DMRgenes<-rbind(DMRgenes, both)
}
write.csv(DMRgenes, 'DMR_H_vs_M_nearest_genes.csv', row.names=FALSE)
write.table(DMRgenes, 'DMR_H_vs_M_nearest_genes.txt', row.names=FALSE)






