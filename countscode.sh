#install HTSeq
#install pip and add to path
#sudo easy_install pip
#sudo nano /etc/paths
#/Users/smartlab/Library/Python/2.7/bin
#install xcode from appstore
#install NumPy
#python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
#install HTSeq
#cd htseq
#python setup.py build
#sudo python setup.py install

#run htseq on smallest bam file as a test
#htseq-count -f bam -r pos -i Parent -o /Volumes/Hyden⁩/blh226⁩/mapped_reads_full_genome/counts/10309.3.156016.ATTACTC-TAATCTT /Volumes/Hyden⁩/blh226⁩/mapped_reads_full_genome/10309.3.156016.ATTACTC-TAATCTT_Aligned.out.sorted.bam /Volumes/Hyden⁩/blh226⁩/Spurpurea⁩/v5.1⁩/annotation⁩/Spurpurea_519_v5.1.gene_exons.gff3
#Previous method outputs unecessary sam file-send output to text file
#htseq-count -f bam -r name -i Parent -m union /Volumes/Hyden⁩/blh226⁩/mapped_reads_full_genome/10309.3.156016.ATTACTC-TAATCTT_Aligned.out.sorted.bam /Volumes/Hyden⁩/blh226⁩/Spurpurea⁩/v5.1⁩/annotation⁩/Spurpurea_519_v5.1.gene_exons.gff3 > /Volumes/Hyden/blh226/mapped_reads_full_genome/counts/10309.3.156016.ATTACTC-TAATCTT_counts_2.txt

#try feature counts


#############################################
#featureCounts
#############################################



#tar zxvf subread-1.6.4-source.tar.gz
#/Volumes/Hyden/blh226/subread-1.6.4-source/src
#make -f Makefile.MacOS
#PATH=$PATH:/Volumes/Hyden/blh226/subread-1.6.4-source/bin
#test file
cd mapped_reads_full_genome
featureCounts -a /Volumes/Hyden/blh226/Spurpurea/v5.1/annotation/Spurpurea_519_v5.1.gene_exons.gff3 -o 10309.3.156016.ATTACTC-TAATCTT_counts.tsv <10309.3.156016.ATTACTC-TAATCTT_Aligned.out.sorted.bam> -F -g 'Parent' -M -p -T 23 --fraction
#worked well, now do for loop to get raw counts for each mapping
#Full Genome
for file in *.sorted.bam
do
	baseFilename=`basename $file _Aligned.out.sorted.bam`
	echo "${baseFilename}"
	featureCounts -a /Volumes/Hyden/blh226/Spurpurea/v5.1/annotation/Spurpurea_519_v5.1.gene_exons.gff3 -o counts/"${baseFilename}"_counts.tsv <"${baseFilename}"_Aligned.out.sorted.bam> -F -g 'Parent' -M -p -T 23 --fraction
	echo 'Thank you, next'
done

#No 15Z
cd /Volumes/Hyden/blh226/mapped_reads_no_15Z
for file in *.sorted.bam
do
	baseFilename=`basename $file _Aligned.out.sorted.bam`
	echo "${baseFilename}"
	featureCounts -a /Volumes/Hyden/blh226/Genome_no_15Z/Spurpurea_519_v5.1.gene_exons_no_15Z.gff3 -o counts/"${baseFilename}"_counts.tsv <"${baseFilename}"_Aligned.out.sorted.bam> -F -g 'Parent' -M -p -T 23 --fraction
	echo 'Thank you, next'
done

#No 15W
cd /Volumes/Hyden/blh226/mapped_reads_no_15W
for file in *.sorted.bam
do
	baseFilename=`basename $file _Aligned.out.sorted.bam`
	echo "${baseFilename}"
	featureCounts -a /Volumes/Hyden/blh226/Genome_no_15W/Spurpurea_519_v5.1.gene_exons_no_15W.gff3 -o counts/"${baseFilename}"_counts.tsv <"${baseFilename}"_Aligned.out.sorted.bam> -F -g 'Parent' -M -p -T 23 --fraction
	echo 'Thank you, next'
done

#Fish Creek with 15W
cd /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W
for file in *.sorted.bam
do
	baseFilename=`basename $file _Aligned.out.sorted.bam`
	echo "${baseFilename}"
	featureCounts -a /Volumes/Hyden/blh226/Genome.FC/SpurpureaFishCreek_518_v3.1.gene_exons_with_15W.gff3 -o counts/"${baseFilename}"_counts.tsv <"${baseFilename}"_Aligned.out.sorted.bam> -F -g 'Parent' -M -p -T 23 --fraction
	echo 'Thank you, next'
done


#########################################
#Combine and clean up files 
####################################
cd /Volumes/Hyden/blh226/mapped_reads_full_genome/counts
R
#read in all files
list<-list.files(pattern="*.tsv")
list2<-c()
for(i in 1:187){
	list2<-c(list2, list[2*i-1])
}
list<-list2 #list is all files minus summary files
allfiles<-lapply(list, read.table)
for(i in 1:length(allfiles)){
	colnames(allfiles[[i]])<-as.character(unlist(allfiles[[i]][1,]))
	allfiles[[i]]<-allfiles[[i]][-1,]
	allfiles[[i]]$Chr<-gsub(";.*","",allfiles[[i]]$Chr)
	allfiles[[i]]$STDIN<-as.numeric(as.character(allfiles[[i]]$STDIN))
} #gives list of all files with first row as column names, removes excess chromosome names

#need to extract only rows that have chr 15Z, 15W
allfiles15Zsubset<-as.character(allfiles[[1]][allfiles[[1]]$Chr=='Chr15Z',1]) #Gives vector with 15Z geneids
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	Zvec<-c()
	for(j in 1:length(temp$Geneid)){
		if(temp[j,2]=='Chr15Z'){
			Zvec<-c(Zvec, temp[j,7]) #Zvec is a vector of 15Z gene counts
		}
	}
	allfiles15Zsubset<-cbind(allfiles15Zsubset, Zvec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfiles15Zsubset)<-c('Geneid', list) #gives appropriate column names
fullgenome15Zcounts<-allfiles15Zsubset
write.csv(fullgenome15Zcounts, 'fullgenome15Zcounts.csv')
#repeat for 15W
allfiles15Wsubset<-as.character(allfiles[[1]][allfiles[[1]]$Chr=='Chr15W',1])
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	Wvec<-c()
	for(j in 1:length(temp$Geneid)){
		if(temp[j,2]=='Chr15W'){
			Wvec<-c(Wvec, temp[j,7]) #Zvec is a vector of 15W gene counts
		}
	}
	allfiles15Wsubset<-cbind(allfiles15Wsubset, Wvec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfiles15Wsubset)<-c('Geneid', list) #gives appropriate column names
fullgenome15Wcounts<-allfiles15Wsubset
write.csv(fullgenome15Wcounts, 'fullgenome15Wcounts.csv')
#repeat for all counts
allfilescounts<-as.character(allfiles[[1]][,1])
chr<-as.character(allfiles[[1]][,2])
allfilescounts<-cbind(allfilescounts, chr)
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	vec<-c()
	for(j in 1:length(temp$Geneid)){
		vec<-c(vec, temp[j,7]) #Zvec is a vector of gene counts
	}
	allfilescounts<-cbind(allfilescounts, vec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfilescounts)<-c('Geneid', 'Chr', list)
fullgenomecounts<-allfilescounts
write.csv(fullgenomecounts, 'fullgenomecounts.csv')
quit()

#repeat process for other alignments
#############################################
	#no 15Z: 15W
cd /Volumes/Hyden/blh226/mapped_reads_no_15Z/counts
R
#read in all files
list<-list.files(pattern="*.tsv")
list2<-c()
for(i in 1:187){
	list2<-c(list2, list[2*i-1])
}
list<-list2 #list is all files minus summary files
allfiles<-lapply(list, read.table)
for(i in 1:length(allfiles)){
	colnames(allfiles[[i]])<-as.character(unlist(allfiles[[i]][1,]))
	allfiles[[i]]<-allfiles[[i]][-1,]
	allfiles[[i]]$Chr<-gsub(";.*","",allfiles[[i]]$Chr)
	allfiles[[i]]$STDIN<-as.numeric(as.character(allfiles[[i]]$STDIN))
} #gives list of all files with first row as column names, removes excess chromosome names

#repeat for 15W
allfiles15Wsubset<-as.character(allfiles[[1]][allfiles[[1]]$Chr=='Chr15W',1])
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	Wvec<-c()
	for(j in 1:length(temp$Geneid)){
		if(temp[j,2]=='Chr15W'){
			Wvec<-c(Wvec, temp[j,7]) #Zvec is a vector of 15W gene counts
		}
	}
	allfiles15Wsubset<-cbind(allfiles15Wsubset, Wvec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfiles15Wsubset)<-c('Geneid', list) #gives appropriate column names
noZ15Wcounts<-allfiles15Wsubset
write.csv(noZ15Wcounts, 'noZ15Wcounts.csv')
#repeat for all counts
allfilescounts<-as.character(allfiles[[1]][,1])
chr<-as.character(allfiles[[1]][,2])
allfilescounts<-cbind(allfilescounts, chr)
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	vec<-c()
	for(j in 1:length(temp$Geneid)){
		vec<-c(vec, temp[j,7]) #Zvec is a vector of gene counts
	}
	allfilescounts<-cbind(allfilescounts, vec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfilescounts)<-c('Geneid', 'Chr', list)
noZcounts<-allfilescounts
write.csv(noZcounts, 'noZcounts.csv')
quit()

##########################################
	#no 15W: 15Z
cd /Volumes/Hyden/blh226/mapped_reads_no_15W/counts
R
#read in all files
list<-list.files(pattern="*.tsv")
list2<-c()
for(i in 1:187){
	list2<-c(list2, list[2*i-1])
}
list<-list2 #list is all files minus summary files
allfiles<-lapply(list, read.table)
for(i in 1:length(allfiles)){
	colnames(allfiles[[i]])<-as.character(unlist(allfiles[[i]][1,]))
	allfiles[[i]]<-allfiles[[i]][-1,]
	allfiles[[i]]$Chr<-gsub(";.*","",allfiles[[i]]$Chr)
	allfiles[[i]]$STDIN<-as.numeric(as.character(allfiles[[i]]$STDIN))
} #gives list of all files with first row as column names, removes excell chromosome names

#repeat for 15Z
allfiles15Zsubset<-as.character(allfiles[[1]][allfiles[[1]]$Chr=='Chr15Z',1])
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	Zvec<-c()
	for(j in 1:length(temp$Geneid)){
		if(temp[j,2]=='Chr15Z'){
			Zvec<-c(Zvec, temp[j,7]) #Zvec is a vector of 15W gene counts
		}
	}
	allfiles15Zsubset<-cbind(allfiles15Zsubset, Zvec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfiles15Zsubset)<-c('Geneid', list) #gives appropriate column names
noW15Zcounts<-allfiles15Zsubset
write.csv(noW15Zcounts, 'noW15Zcounts.csv')
#repeat for all counts
allfilescounts<-as.character(allfiles[[1]][,1])
chr<-as.character(allfiles[[1]][,2])
allfilescounts<-cbind(allfilescounts, chr)
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	vec<-c()
	for(j in 1:length(temp$Geneid)){
		vec<-c(vec, temp[j,7]) #Zvec is a vector of gene counts
	}
	allfilescounts<-cbind(allfilescounts, vec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfilescounts)<-c('Geneid', 'Chr', list)
noWcounts<-allfilescounts
write.csv(noWcounts, 'noWcounts.csv')
quit()

################################################################
	#FC: 15W
cd /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/counts
R
#read in all files
list<-list.files(pattern="*.tsv")
list2<-c()
for(i in 1:187){
	list2<-c(list2, list[2*i-1])
}
list<-list2 #list is all files minus summary files
allfiles<-lapply(list, read.table)
for(i in 1:length(allfiles)){
	colnames(allfiles[[i]])<-as.character(unlist(allfiles[[i]][1,]))
	allfiles[[i]]<-allfiles[[i]][-1,]
	allfiles[[i]]$Chr<-gsub(";.*","",allfiles[[i]]$Chr)
	allfiles[[i]]$STDIN<-as.numeric(as.character(allfiles[[i]]$STDIN))
} #gives list of all files with first row as column names, removes excell chromosome names
#repeat for 15W
allfiles15Wsubset<-as.character(allfiles[[1]][allfiles[[1]]$Chr=='Chr15W',1])
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	Wvec<-c()
	for(j in 1:length(temp$Geneid)){
		if(temp[j,2]=='Chr15W'){
			Wvec<-c(Wvec, temp[j,7]) #Zvec is a vector of 15W gene counts
		}
	}
	allfiles15Wsubset<-cbind(allfiles15Wsubset, Wvec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfiles15Wsubset)<-c('Geneid', list) #gives appropriate column names
FC15Wonly<-allfiles15Wsubset
write.csv(FC15Wonly, 'FC15Wonly.csv')
#repeat for all counts
allfilescounts<-as.character(allfiles[[1]][,1])
chr<-as.character(allfiles[[1]][,2])
allfilescounts<-cbind(allfilescounts, chr)
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	vec<-c()
	for(j in 1:length(temp$Geneid)){
		vec<-c(vec, temp[j,7]) #Zvec is a vector of gene counts
	}
	allfilescounts<-cbind(allfilescounts, vec) #allfiles15Zsubset is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
colnames(allfilescounts)<-c('Geneid', 'Chr', list)
FCand15Wcounts<-allfilescounts
write.csv(FCand15Wcounts, 'FCand15Wcounts.csv')
quit()
	#should have five total files (three 15W and two 15Z)

#generate raw counts and RPKM 
	#install and test cufflinks 
#tar zxvf /Volumes/Hyden/blh226/boost_1_70_0.tar.gz 
#cd boost_1_70_0
#bash bootstrap.sh
#./b2
#./bjam --toolset=darwin architecture=x86 link=static runtime-link=static --layout=versioned stage install
#cd .. 
#tar zxvf eigen-eigen-323c052e1731.tar.gz
#cd eigen-eigen-323c052e1731
#mv /Volumes/Hyden/blh226/eigen-eigen-323c052e1731 /usr/local/include
#tar zxvf /Volumes/Hyden/blh226/cufflinks-2.1.1.OSX_x86_64.tar.gz
#export PATH=$PATH:/Volumes/Hyden/blh226/cufflinks-2.1.1.OSX_x86_64
#PATH=$PATH:/Volumes/Hyden/blh226/boost_1_70_0
#cd ..
#cufflinks ./test_data.sam
	#test cufflinks on a single dataset
#cufflinks -o /Volumes/Hyden/mapped_reads_final/RPKM -p 24 --library-type fr-firststrand -m 292 -v -g Genome_final/Spurpurea_519_v5.1.gene_exons_187Zspecific.gff3 mapped_reads_final/mapped_reads_final_merged_subset/LIB-10X-317-002-1_Aligned.out.sorted.bam
#cufflinks slow and didn't work 



##################################################
#featureCounts 
################################################
featureCounts -a /Volumes/Hyden/blh226/Genome_final/Spurpurea_519_v5.1.gene_exons_187Zspecific.gff3 -o /Volumes/Hyden/blh226/mapped_reads_final/counts/10X-317-002_counts.tsv </Volumes/Hyden/blh226/mapped_reads_final/mapped_reads_final_merged_subset/1LIB-10X-317-002-1_Aligned.out.sorted.bam> -F -g Parent -M -p -T 24 --fraction

cd mapped_reads_final/mapped_reads_final_merged_subset
for file in *.sorted.bam
do
	echo $file
	baseFilename=`basename $file Aligned.out.sorted.bam`
	echo "${baseFilename}"
	featureCounts -a /Volumes/Hyden/blh226/Genome_final/Spurpurea_519_v5.1.gene_exons_187Zspecific.gff3 -o /Volumes/Hyden/blh226/mapped_reads_final/counts/"${baseFilename}"counts.tsv </Volumes/Hyden/blh226/mapped_reads_final/mapped_reads_final_merged_subset/"${baseFilename}"Aligned.out.sorted.bam> -F -g Parent -M -p -T 24 --fraction
	echo 'Thank you, next'
done

##############################################
#Clean up tsv files and combine data 
###############################################
cd /Volumes/Hyden/blh226/mapped_reads_final/counts
R
#read in all files
list<-list.files(pattern="*.tsv")
list2<-c()
for(i in 1:183){
	list2<-c(list2, list[2*i-1])
}
list<-list2 #list is all files minus summary files
allfiles<-lapply(list, read.table)
for(i in 1:length(allfiles)){
	colnames(allfiles[[i]])<-as.character(unlist(allfiles[[i]][1,]))
	allfiles[[i]]<-allfiles[[i]][-1,]
	allfiles[[i]]$Chr<-gsub(";.*","",allfiles[[i]]$Chr)
	allfiles[[i]]$STDIN<-as.numeric(as.character(allfiles[[i]]$STDIN))
} #gives list of all files with first row as column names, removes excess chromosome names

#combine all files into single table (note, nested for loop takes a long time)
allfilescombined<-as.character(allfiles[[1]]$Chr) #Gives vector with geneids
for(i in 1:length(allfiles)){
	temp<-allfiles[[i]] #temp has one file in it
	vec<-c() #holding variable 
	for(j in 1:length(temp$Geneid)){
		vec<-c(vec, temp[j,7]) #vec is a vector of gene counts
	}
	allfilescombined<-cbind(allfilescombined, vec) #allfilescombined is a dataframe containing geneid in the first column, and counts by sample in the remaining columns
}
#create header with sample/file names 
colnames(allfilescombined)<-c('Geneid', list) #gives appropriate column names
write.csv(allfilescombined, 'final_counts_combined.csv')

