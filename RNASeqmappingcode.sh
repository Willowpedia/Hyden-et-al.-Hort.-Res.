cd /Volumes/Hyden/blh226
#Download binaries for bowtie2
#wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip
#unzip bowtie2-2.2.4-linux-x86_64.zip
#download binaries for tophat2
#wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.13.Linux_x86_64.tar.gz
#tar -xzf tophat-2.0.13.Linux_x86_64.tar.gz

#index full genome
#bowtie2-build /Volumes/Hyden/blh226/Spurpurea/v5.1/assembly/Spurpurea_519_v5.0.fa Spurpurea_519_v5.0
#decouple fastq file
#deinterleave_fastq.sh < /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG.anqrpht.fastq /Volumes/Hyeden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq
#test alignment of one RNASeq fastq file
#tophat -r 200 Spurpurea_519_v5.0 /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq
#took over ten hours to map, but results look pretty good
###############################################
#need to specify options
#check that names of index match GFF3 file
#bowtie2-inspect --names Spurpurea_519_v5.0
#run tophat
#tophat -o /Volumes/Hyden/blh226/test_output -r 200 -p 24 -G /Volumes/Hyden/blh226/Spurpurea/v5.1/annotation/Spurpurea_519_v5.1.gene_exons.gff3 Spurpurea_519_v5.0 /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq
#failed after one hour, error reading SAM file, retried, same error :
#segment-based junction search failed with err =1
#Error opening SAM file /Volumes/HAL9000/blh226/test_output/tmp/left_kept_reads.m2g_um_seg3.bam
#tophat -r 200 -p 24  Spurpurea_519_v5.0 /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq
#also failed, same error

############################################
#Try using bbmap
#cd ..
#install
#tar -xvzf BBMap_38.50b.tar.gz
#test installation
#bbmap/stats.sh in=bbmap/resources/phix174_ill.ref.fa.gz
#index genome
#bbmap ref=/Volumes/Hyden/blh226/Spurpurea/V5.1/assembly/Spurpurea_519_v5.0.fa
#RNAseq mapping basic
#bbmap interleaved=true in=fastq/10262.1.154566.ATTACTC-ATAGAGG.anqrpht.fastq out=10262.1.154566.ATTACTC-ATAGAGG.anqrpht.sam
#RNAseq mapping further specifications
#bbmap -Xmx54g interleaved=true in=fastq/10262.1.154566.ATTACTC-ATAGAGG.anqrpht.fastq out=10262.1.154566.ATTACTC-ATAGAGG.anqrpht.sam averagepairdist=200 rpkm=10262.1.154566.ATTACTC-ATAGAGG.anqrpht.txt threads=22
#took too long
##############################################
bash deinterleave_fastq.sh < /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG.anqrpht.fastq /Volumes/Hyeden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq
#maybe also try STAR
#install
#git clone https://github.com/alexdobin/STAR.git
#cd STAR/source
#run STAR
#index genome only once
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --genomeFastaFiles /Volumes/Hyden/blh226/Spurpurea/V5.1/assembly/Spurpurea_519_v5.0.fa --sjdbGTFfile /Volumes/Hyden/blh226/Spurpurea/v5.1/annotation/Spurpurea_519_v5.1.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 150
#mapping
STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq --outSAMtype BAM Unsorted SortedByCoordinate
#failed not enough disk space
#try without sorting BAM
STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/10262.1.154566.ATTACTC-ATAGAGG_2.anqrpht.fastq --outSAMtype BAM Unsorted
###################################################
#to view bam file (optional)
samtools view -h /Volumes/Hyden/blh226/Aligned.out.bam> Aligned.out.sam
#Need to sort and index bam file
samtools sort -T /tmp/Aligned.out.sorted -o Aligned.out.sorted.bam -@ 22 Aligned.out.bam
samtools index Aligned.out.sorted.bam
#gzip only uses one core, to use all cores use pigz
brew install pigz

#optimize for doing all samples
#first need to deinterleave fastq files-keep zipped to speed up run time
reformat.sh -Xmx56g in=/Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA.anqrpht.fastq.gz out1=/Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_1.anqrpht.fastq.gz out2=/Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_2.anqrpht.fastq.gz
#unzip deinterleaved files, use pigz to increase time (uses all available cores by default)
pigz -d /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_2.anqrpht.fastq.gz
#run STAR for RNAseq mapping using decompressed files, for some reason, trying to sort BAM using STAR doesn't work
STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/11801.1.220153.ATAGCGG-ACCGCTA_ --outReadsUnmapped Fasta --outSAMtype BAM Unsorted
#need to seperately sort BAM file
samtools sort -T /tmp/11801.1.220153.ATAGCGG-ACCGCTA_Aligned.out.sorted -o 11801.1.220153.ATAGCGG-ACCGCTA_Aligned.out.sorted.bam -@ 22 11801.1.220153.ATAGCGG-ACCGCTA_Aligned.out.bam
#Need to index BAM file
samtools index 11801.1.220153.ATAGCGG-ACCGCTA_Aligned.out.sorted.bam
#re-compress deinterleaved fastq files to save storage space
pigz /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/11801.1.220153.ATAGCGG-ACCGCTA_2.anqrpht.fastq
#tell me when its done
echo "11801.1.220153.ATAGCGG-ACCGCTA.anqrpht.fastq finished"

###################################################################################################
#deinterleave 
##################################################################################################

#command to de interleave all files and place in new directory, only need to do once for all subsequent mapping attempts
for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	reformat.sh -Xmx56g in=$file out1=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz out2=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
done

#Repeat for files with .fq.gz extension
for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	reformat.sh -Xmx56g in=$file out1=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz out2=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
done
#repeat for files with .catkin.fastq.gz extension
for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	reformat.sh -Xmx56g in=$file out1=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz out2=/Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
done

##################################
#Mapping
###################################


#command to decompress files map reads, and recompress
for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_ --outReadsUnmapped Fasta --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq
	echo finished "${baseFilename}"
done

for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_ --outReadsUnmapped Fasta --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq
done

for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/STAR_genome_directory --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_ --outReadsUnmapped Fasta --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_full_genome/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done

#Repeat indexing and mapping for genome without 15Z
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir /Volumes/Hyden/blh226/Genome_no_15Z --genomeFastaFiles /Volumes/Hyden/blh226/Genome_no_15Z/Spurpurea_519_v5.0_no_15Z.fa --sjdbGTFfile /Volumes/Hyden/blh226/Genome_no_15Z/Spurpurea_519_v5.1.gene_exons_no_15Z.gff3 --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 150

for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15Z --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq
done

for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15Z --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq
done

for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15Z --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15Z/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done

#Repeat indexing and mapping for genome without 15W
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir /Volumes/Hyden/blh226/Genome_no_15W --genomeFastaFiles /Volumes/Hyden/blh226/Genome_no_15W/Spurpurea_519_v5.0_no_15W.fa --sjdbGTFfile /Volumes/Hyden/blh226/Genome_no_15W/Spurpurea_519_v5.1.gene_exons_no_15W.gff3 --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 150

for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq
done

for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done

for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_no_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_no_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq
done

#Repeat indexing and mapping for Fish Creek with 15W (need to specify extra chromosome seperately)
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir /Volumes/Hyden/blh226/Genome_FishCreek_with_15W --genomeFastaFiles /Volumes/Hyden/blh226/Genome.FC/Salix_purpurea_var_Fishcreek.mainGenome.fasta /Volumes/Hyden/blh226/Genome.FC/Chr15W.fasta --sjdbGTFfile /Volumes/Hyden/blh226/Genome.FC/SpurpureaFishCreek_518_v3.1.gene_exons_with_15W.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150

for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_FishCreek_with_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq
done

for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_FishCreek_with_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done

for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
	STAR --runThreadN 22 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_FishCreek_with_15W --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam -@ 22 /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_FishCreek_with_15W/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq
done

#####################################
##########final mapping##############
####################################
#index genome only once
cd /Volumes/Hyden/blh226
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /Volumes/Hyden/blh226/Genome_final --genomeFastaFiles /Volumes/Hyden/blh226/Genome_final/Spurpurea_519_v5.0_187Zspecific.fa --sjdbGTFfile /Volumes/Hyden/blh226/Genome_final/Spurpurea_519_v5.1.gene_exons_187Zspecific.gff3 --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 150

#mapping 

for file in fastq/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 24 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_final --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_ --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam -@ 24 /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done

#will need to re-do 10262.1.154566.CGCTCAT-AGGCGAA, de-interleaved file was truncated for some reason,need to recreate file and re-run mapping
pigz -d fastq/10262.1.154566.CGCTCAT-AGGCGAA.anqrpht.fastq

reformat.sh -Xmx56g in=/Volumes/Hyden/blh226/fastq/10262.1.154566.CGCTCAT-AGGCGAA.anqrpht.fastq out1=/Volumes/Hyden/blh226/fastq/10262.1.154566.CGCTCAT-AGGCGAA_1.anqrpht.fastq.gz out2=/Volumes/Hyden/blh226/fastq/10262.1.154566.CGCTCAT-AGGCGAA_2.anqrpht.fastq.gz

pigz -d fastq/10262.1.154566.CGCTCAT-AGGCGAA.anqrpht.fastq 

pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_2.anqrpht.fastq.gz
STAR --runThreadN 24 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_final --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_final/10262.1.154566.CGCTCAT-AGGCGAA_ --outSAMtype BAM Unsorted
samtools sort -T /tmp/10262.1.154566.CGCTCAT-AGGCGAA_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_final/10262.1.154566.CGCTCAT-AGGCGAA_Aligned.out.sorted.bam -@ 24 /Volumes/Hyden/blh226/mapped_reads_final/10262.1.154566.CGCTCAT-AGGCGAA_Aligned.out.bam
samtools index /Volumes/Hyden/blh226/mapped_reads_final/10262.1.154566.CGCTCAT-AGGCGAA_Aligned.out.sorted.bam
pigz /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/10262.1.154566.CGCTCAT-AGGCGAA_2.anqrpht.fastq

#also need to do duplicate files 
for file in fastq/duplicates/*.anqrpht.fastq.gz
do
	baseFilename=`basename $file .anqrpht.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq.gz
	STAR --runThreadN 24 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_final --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_ --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam -@ 24 /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.anqrpht.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.anqrpht.fastq
done
#fq.gz
for file in fastq/*.fq.gz
do
	baseFilename=`basename $file .fq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq.gz
	STAR --runThreadN 24 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_final --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_ --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam -@ 24 /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.fq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.fq
done
#catkin.fastq.gz 
for file in fastq/*.catkin.fastq.gz
do
	baseFilename=`basename $file .catkin.fastq.gz`
	echo "${baseFilename}"
	pigz -d /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq.gz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq.gz
	STAR --runThreadN 24 --runMode alignReads --genomeDir /Volumes/Hyden/blh226/Genome_final --readFilesIn /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq --outFileNamePrefix /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_ --outSAMtype BAM Unsorted
	samtools sort -T /tmp/"${baseFilename}"_Aligned.out.sorted -o /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam -@ 24 /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.bam
	samtools index /Volumes/Hyden/blh226/mapped_reads_final/"${baseFilename}"_Aligned.out.sorted.bam
	pigz /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_1.catkin.fastq /Volumes/Hyden/blh226/fastq/deinterleaved/"${baseFilename}"_2.catkin.fastq
done


######################################################################################################################################################
#generate some mapping statistics, including comparing between methods to make sure everything is working prior to moving forward.  
#######################################################################################################################################################


for file in mapped_reads_full_genome/*.final.out
do 
	sed -n 10p $file | grep -o -E '[0-9.]+' >> mapped_reads_full_genome/unique_mapped_full.txt
done 

for file in mapped_reads_no_15Z/*.final.out
do 
	sed -n 10p $file | grep -o -E '[0-9.]+' >> mapped_reads_no_15Z/unique_mapped_no15Z.txt
done 

for file in mapped_reads_no_15W/*.final.out
do 
	sed -n 10p $file | grep -o -E '[0-9.]+' >> mapped_reads_no_15W/unique_mapped_no15W.txt
done 

for file in mapped_reads_final/*.final.out
do 
	sed -n 10p $file | grep -o -E '[0-9.]+' >> mapped_reads_final/unique_mapped_final.txt
done 

R 
setwd('/Volumes/Hyden/blh226')
full<-read.table('mapped_reads_full_genome/unique_mapped_full.txt')
noZ<-read.table('mapped_reads_no_15Z/unique_mapped_no15Z.txt')
noW<-read.table('mapped_reads_no_15W/unique_mapped_no15W.txt')
final<-read.table('mapped_reads_final/unique_mapped_final.txt')
full<-unlist(full)
noW<-unlist(noW)
noZ<-unlist(noZ)
final<-unlist(final)
mean(full)
sd(full)
mean(noZ)
sd(noZ)
mean(noW)
sd(noW)
mean(final)
sd(final)
hist(full, main='full genome', col='gold', xlab='Percent Unique Mapped Reads')
hist(noZ, main='full genome minus 15Z', col='gold', xlab='Percent Unique Mapped Reads')
hist(noW, main='full genome minus 15W', col='gold', xlab='Percent Unique Mapped Reads')
hist(final, main='full genome with Z scaffold (final)', col='gold', xlab='Percent Unique Mapped Reads')


#########################################################
#rename files according to sample name 
#########################################################
cd mapped_reads_final
for file in *_Log.final.out
do
	baseFilename=`basename $file _Log.final.out`
	echo "${baseFilename}"
	awk '$1 == "'${baseFilename}'" {print $3}' /Volumes/Hyden/blh226/fastq/female_male_hermie2.keyfile
	newname=`awk '$1 == "'${baseFilename}'" {print $3}' /Volumes/Hyden/blh226/fastq/female_male_hermie2.keyfile`
	echo "${newname}"
	rename "s/$baseFilename/$newname"/ *
done
#samtools merge to merge bams from duplicate files 
for file in *-2_Aligned.out.sorted.bam 
do
	baseFilename=`basename $file -2_Aligned.out.sorted.bam`
	echo "${baseFilename}"
	samtools merge "${baseFilename}"_merged_samples.Aligned.out.sorted.bam "${baseFilename}"-2_Aligned.out.sorted.bam "${baseFilename}"-1_Aligned.out.sorted.bam -@ 24
	samtools index "${baseFilename}"_merged_samples.Aligned.out.sorted.bam
done

samtools merge LIB-94003_merged_samples.Aligned.out.sorted.bam LIB-00-22-002-1_Aligned.out.sorted.bam LIB-94003-1_Aligned.out.sorted.bam -@ 24
samtools index LIB-94003_merged_samples.Aligned.out.sorted.bam
samtools merge LIB-9882-34_merged_samples.Aligned.out.sorted.bam LIB-9882-34_R1-1_Aligned.out.sorted.bam LIB-9882-34_R2-1_Aligned.out.sorted.bam LIB-9882-34_R3-1_Aligned.out.sorted.bam -@ 24
samtools index LIB-9882-34_merged_samples.Aligned.out.sorted.bam
#manually move final bam files into new directory for further analysis 
