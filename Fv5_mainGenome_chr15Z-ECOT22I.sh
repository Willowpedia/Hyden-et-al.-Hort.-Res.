#########################################################################################################################
#########################################################################################################################
#
# Project: Willow GBS Pipeline script
# Reference genome: Female 'Salix purpurea' clone 94006 v4 (PacBio polished assembly)
# Population: F2 Salix purpurea family 317 (94006 x 94001 -> F1 family (9882-41 x 9882-34) -> F2 family)
# Restriction enzyme: EcoT22I
# Date: 2018-03-16
#
# Created by:
#
# Craig H. Carlson
# 
# Horticulture Section
# School of Integrative Plant Sciences
# Cornell University
# NYS Agricultural Experiment Station
# 630 West North Street
# Geneva, NY 14456
#
#########################################################################################################################
#################################################### Begin GBS Pipeline #################################################
# Normalize fasta
picard NormalizeFasta I=Reference_genome/Fv5_main_genome_chr15Z.fasta O=Reference_genome/Fv5_mainGenome_chr15Z.fasta

# Build a BWA-index for the reference genome
bwa index -a is Reference_genome/Fv5_mainGenome_chr15Z.fasta

# Build a SAMTools index for the reference genome
samtools faidx Reference_genome/Fv5_mainGenome_chr15Z.fasta
#cd "Hyden/F2_EcoT22I"

# FastqToTagCountPlugin 
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -FastqToTagCountPlugin -i fastq -k AoB.EcoT22I.keyfile -e EcoT22I -s 350000000 -c 1 -o tagCounts -endPlugin -runfork1 > tagCounts/FastqToTagCountPlugin.log 2> tagCounts/FastqToTagCountPlugin.err

# MergeMultipleTagCountPlugin 
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeMultipleTagCountPlugin -i tagCounts -o mergedTagCounts/mergedTagCounts.cnt -c 5 -endPlugin -runfork1 > mergedTagCounts/mergedTagCounts_c5_nb.log 2> mergedTagCounts/mergedTagCounts_c5_nb.err

# BinaryToTextPlugin
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -BinaryToTextPlugin -i mergedTagCounts/mergedTagCounts.cnt -o mergedTagCounts/mergedTagCounts.cnt.txt -t tagCounts -endPlugin -runfork1

# TagCountToFastqPlugin
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -TagCountToFastqPlugin -i mergedTagCounts/mergedTagCounts.cnt -o mergedTagCounts/mergedTagCounts.fq -c 5 -endPlugin -runfork1 > mergedTagCounts/mergedTagCounts_c5_cvtnb.log 2> mergedTagCounts/mergedTagCounts_c5_cvtnb.err

# FastqToTBTPlugin
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -FastqToTBTPlugin -i fastq -k AoB.EcoT22I.keyfile -e EcoT22I -o tbt -y -t mergedTagCounts/mergedTagCounts.cnt -endPlugin -runfork1 > tbt/fq2tbt.log 2> tbt/fq2tbt.err

# MergeTagsByTaxaFilesPlugin - No merged Taxa
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeTagsByTaxaFilesPlugin -i tbt -s 350000000 -o mergedTBT/mergedTBT.tbt.byte -endPlugin -runfork1 > mergedTBT/tbt_nb.log 2> mergedTBT/tbt_nb.err

# MergeTagsByTaxaFilesPlugin - Merged Taxa
# Use "-x" to merge tag counts of taxa with identical short names
#~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeTagsByTaxaFilesPlugin -i tbt -x -s 350000000 -o mergedTBT/mergedTBT.x.tbt.byte -endPlugin -runfork1 > mergedTBT/tbt_x_nb.log 2> mergedTBT/tbt_x_nb.err

####################################################
# Map mergedTagCounts to reference (-o SAM)
####################################################

# BWA-mem Perl Script (-i Fastq [.fq] -o SequenceAlignment [.sam])
bwa mem -t 12 Reference_genome/Fv5_mainGenome_chr15Z.fasta mergedTagCounts/mergedTagCounts.fq > topm/Fv5_mainGenome_chr15Z.eco.sam

# First run CHMOD on that Perl shiznit ("chmod +x" allows us to execute the script)!
chmod +x filter_bwa_mem.pl

# Convert from BWA-format to a readable Bowtie2-format SAM file
perl ./filter_bwa_mem.pl topm/Fv5_mainGenome_chr15Z.eco.sam topm/Fv5_mainGenome_chr15Z.eco.filt.sam 1

# SAMConverterPlugin - BWA
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -SAMConverterPlugin -i topm/Fv5_mainGenome_chr15Z.eco.filt.sam -o topm/Fv5_mainGenome_chr15Z.eco.topm -endPlugin -runfork1 > topm/topm_Z.log 2> topm/topm_Z.err

####################################################
# "Taxa UNMERGED" - Call variants (-o VCF) 
####################################################

# tbt2vcfPlugin: "UNMERGED" TBT (-i mergedTBT/mergedTBT.tbt.byte)
~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -tbt2vcfPlugin -i mergedTBT/mergedTBT.tbt.byte -m topm/Fv5_mainGenome_chr15Z.eco.topm -o tbt2vcf15Z -mnMAF 0.05 -mnLCov 0.0 -ak 4 -s 1 -e 844 -endPlugin -runfork1 > tbt2vcf15Z/tbt2vcf.log 2> tbt2vcf15Z/tbt2vcf.err

# Set working directory to tbt2vcf output
#cd "Hyden/F2_EcoT22I"

# Compress the mergedTBT for each chromosome
gzip tbt2vcf15Z/mergedTBT.c*

# Give vcftools the path to PERL5
export PERL5LIB=~/programs/vcftools_0.1.13/lib/perl5/site_perl

# Concatenate VCFs and gzip the output
# ~/programs/vcftools_0.1.13/bin/vcf-concat tbt2vcf/mergedTBT.c*.gz | gzip -c > vcf/Fv5_mainGenome_chr15Z.eco.vcf.gz
~/programs/vcftools_0.1.13/bin/vcf-concat tbt2vcf15Z/mergedTBT.c*.gz | gzip -c > hyden_vcf/Fv5_mainGenome_chr15Z.eco.vcf.gz

# De-compress VCF to merge and sort
# gunzip -cd vcf/Fv5_mainGenome_chr15Z.eco.vcf.gz > vcf/Fv5_mainGenome_chr15Z.eco.vcf
gunzip -cd hyden_vcf/Fv5_mainGenome_chr15Z.eco.vcf.gz > hyden_vcf/Fv5_mainGenome_chr15Z.eco.vcf

# Merge duplicate SNPs
# ~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -MergeDuplicateSNP_vcf_Plugin -i vcf/Fv5_mainGenome_chr15Z.eco.vcf -ak 4 -o vcf/Fv5_mainGenome_chr15Z.eco.mergedDup.vcf -endPlugin -runfork1
~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -MergeDuplicateSNP_vcf_Plugin -i hyden_vcf/Fv5_mainGenome_chr15Z.eco.vcf -ak 4 -o hyden_vcf/Fv5_mainGenome_chr15Z.eco.mergedDup.vcf -endPlugin -runfork1

# Sort the VCF
# ~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -SortGenotypeFilePlugin -inputFile vcf/AFv5_mainGenome_chr15Z.eco.mergedDup.vcf -outputFile vcf/Fv5_mainGenome_chr15Z.eco.mergedDupSorted.vcf -fileType VCF -endPlugin -runfork1
~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -SortGenotypeFilePlugin -inputFile hyden_vcf/Fv5_mainGenome_chr15Z.eco.mergedDup.vcf -outputFile hyden_vcf/Fv5_mainGenome_chr15Z.eco.mergedDupSorted.vcf -fileType VCF -endPlugin -runfork1

# Compress and index sorted VCF
bgzip hyden_vcf/Fv5_mainGenome_chr15Z.eco.MergedDupSorted.vcf
tabix -p vcf hyden_vcf/Fv5_mainGenome_chr15Z.eco.MergedDupSorted.vcf.gz

# TASSEL GUI
~/programs/tassel-5-standalone/start_tassel.pl -Xms2g -Xmx56g

# Compress and index sorted VCFs
# for VCF in $(ls *.vcf); do bgzip $VCF; tabix -p vcf $VCF.gz; done

####################################################
# "Taxa MERGED" - Call variants (-o VCF) 
####################################################

# tbt2vcfPlugin: "MERGED" TBTx (-i mergedTBT/mergedTBTx.tbt.byte)
#~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx64g -fork1 -tbt2vcfPlugin -i mergedTBT/mergedTBT.x.tbt.byte -m topm/Salix_purpurea_var_94006_v4.mainGenome.eco.topm -o tbtx2vcf -mnMAF 0.1 -mnLCov 0.2 -ak 4 -s 1 -e 363 -endPlugin -runfork1 > tbtx2vcf/tbtx2vcf.log 2> tbtx2vcf/tbtx2vcf.err

# Set working directory to tbt2vcf output
#cd "/Volumes/HAL9000/genomics/salix_ngs_fastq/fastq_gbs/F2_EcoT22I"

# Compress the mergedTBT for each chromosome
#gzip tbtx2vcf/mergedTBT.c*

# Give vcftools the path to PERL v5
#export PERL5LIB=~/programs/vcftools_0.1.13/lib/perl5/site_perl

# Concatenate VCFs and gzip the output
#~/programs/vcftools_0.1.13/bin/vcf-concat tbtx2vcf/mergedTBT.c*.gz | gzip -c > vcf_mergedTaxa/AoB_94006_v4.mainGenome.eco.vcf.gz

# De-compress VCF to merge and sort
#gunzip -cd vcf_mergedTaxa/AoB_94006_v4.mainGenome.eco.vcf.gz > vcf_mergedTaxa/AoB_94006_v4.mainGenome.eco.vcf

# Merge duplicate SNPs
#~/programs/tassel3-standalone/run_pipeline.pl -Xmx56g -fork1 -MergeDuplicateSNP_vcf_Plugin -i vcf_mergedTaxa/AoB_94006_v4.mainGenome.eco.vcf -ak 4 -o vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDup.eco.vcf -endPlugin -runfork1

# Sort the VCF
#~/programs/tassel-5-standalone/run_pipeline.pl -Xmx56g -fork1 -SortGenotypeFilePlugin -inputFile vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDup.eco.vcf -outputFile vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDupSorted.eco.vcf -fileType VCF -endPlugin -runfork1

#rm vcf/AoB_94006_v4.mainGenome.eco.vcf
#rm vcf/AoB_94006_v4.mainGenome.eco.mergedDup.vcf

# awk '$0="eco."$0' vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDupSorted.eco.vcf > vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDupSorted.ecoNames.vcf

# Compress and index sorted VCF
#bgzip vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDupSorted.eco.vcf
#tabix -p vcf vcf_mergedTaxa/AoB_94006_v4.mainGenome.mergedDupSorted.eco.vcf.gz

# TASSEL GUI
#~/programs/tassel-5-standalone/start_tassel.pl -Xms2g -Xmx56g

################################################## End GBS Pipeline #####################################################
#########################################################################################################################

# Validate VCF
# vcf-validator vcf/AoB_94006_v4.mainGenome.mergedDupSorted.vcf.gz

####################################################
# Annotate VCF with GZIP
####################################################

# CHR START END   ANNOTATION 
# 1   12345 22345 gene1 
# 1   67890 77890 gene2 

#bgzip Salix_purpurea_var_94006_v4_gene_annotations.gff
#tabix -s 1 -b 2 -e 3 Salix_purpurea_var_94006_v4_gene_annotations.gff.gz

#bcftools annotate -a Salix_purpurea_var_94006_v4_gene_annotations.gff.gz -c CHR,START,END,GENE -O z -h key=INFO,ID=ANN,Number=1,Type=Integer,Description='Salix purpurea var 94006 v4 gene annotation' --threads 12 -o common_parent_gbs_94006v4_mergedDupSorted.gff.vcf common_parent_gbs_94006v4_mergedDupSorted.vcf

####################################################
# Accessory Hapmap functions                     
####################################################

# TagsToSNPByAlignmentPlugin (-i [.tbt.byte] -o [.hmp] -m [.topm] -mUpd [withVariants.topm])
#~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -TagsToSNPByAlignmentPlugin -i mergedTBT/mergedTBT.x.tbt.byte -ref genome/Salix_purpurea_var_94006_v4.mainGenome.fasta -m topm/Salix_purpurea_var_94006_v4.mainGenome.eco.topm -y -mnMAF 0.1 -sC 1 -eC 363 -o hapmap/AoB_94006_v4.mainGenome.eco.chr+.hmp.txt -endPlugin -runfork1

# MergeIdenticalTaxaPlugin APEKI
~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx12g -fork1 -h AoB_94006_v4.mainGenome.ape.filtTaxa.hmp.txt -separate -export -runfork1

~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -MergeIdenticalTaxaPlugin -hmp AoB_94006_v4.mainGenome.ape_chrom+.hmp.txt -o AoB_94006_v4.mainGenome.ape_mergedTaxa_chrom+.hmp.txt -hetFreq 0.75 -sC 1 -eC 355 -endPlugin -runfork1

# MergeIdenticalTaxaPlugin ECOT22I
#~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx42g -fork1 -h AoB_94006_v4.mainGenome.ape.eco.filtTaxa.hmp.txt -separate -export -runfork1

#~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -MergeIdenticalTaxaPlugin -hmp AoB_94006_v4.mainGenome.ape.eco_chrom+.hmp.txt -o AoB_94006_v4.mainGenome.ape.eco_mergedTaxa_chrom+.hmp.txt -hetFreq 0.75 -sC 1 -eC 355 -endPlugin -runfork1

#dat <- fread("AoB_94006_v4.mainGenome.ape.eco.filtTaxa.hmp.txt", sep="\t", na.strings=NULL, check.names=F, data.table=F)
#datNames <- names(dat)
#newDatNames <- as.character(gsub("_", "-", datNames))
#names(dat) <- newDatNames

#write.table(dat, "AoB_94006_v4.mainGenome.ape.eco.unmergedTaxa.hmp.txt", sep="\t", row.names=F, quote=F)

#~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx42g -fork1 -h AoB_94006_v4.mainGenome.ape.eco.unmergedTaxa.syn.hmp.txt -separate -export -runfork1

#~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -MergeIdenticalTaxaPlugin -hmp AoB_94006_v4.mainGenome.ape.eco.unmergedTaxa_chrom+.hmp.txt -o AoB_94006_v4.mainGenome.ape.eco_mergedTaxa_chrom+.hmp.txt -hetFreq 0.75 -sC 1 -eC 355 -endPlugin -runfork1

# MergeDuplicateSNPsPlugin (-i HapMap [.hmp] -o HapMap_MergedSNPs [.hmp])
# ~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -MergeDuplicateSNPsPlugin -hmp hapmap/Spur_94006_v3_taxaMerged_chr+.hmp.txt -callHets -o hapmap/Spur_94006_v3_taxaSNPsMerged_chr+.hmp.txt –callHets -sC 1 -eC 246 -endPlugin -runfork1

# GBSHapMapFiltersPlugin
# ~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -GBSHapMapFiltersPlugin -hmp hapmap/mergedSNPs/Spur_94006_v3_taxaSNPsMerged_chr+.hmp.txt -o hapmap/filt/Spur_94006_v3_taxaSNPsMergedFilt_chr+.hmp.txt -mnTCov 0.01 -mnSCov 0.2 -mnMAF 0.01 -hLD -mnR2 0.2 -mnBonP 0.005 -sC 1 -eC 363 -endPlugin -runfork1

#########################################################################################################################

####################################################
# Plugin Descriptions (in order of appearance)   
####################################################

# 1 - FastqToTagCountPlugin

# Description:  Derives a tagCount list for each FASTQ file in the input directory. Keeps only good reads having a barcode and a cut site and no N's in the useful part of the sequence. Trims off the barcodes and truncates sequences that (1) have a second cut site, or (2) read into the common adapter.
# Input: (1) Barcode key file and (2) FASTQ directory.
# Output: (1) Directory containing a tagCount (.cnt) file for every FASTQ file in the input directory.

# 2 - MergeMultipleTagCountPlugin

# Description:  Merges each tagCount file in the input directory into a single “master” tagCount list. Only keeps tags with a total count (after merger) greater than or equal to that specified by the -c option (minimum number of times a tag must be present to be output).
# Input: Input directory (folder) containing tagCount (.cnt) files.
# Output: Merged tagCount file.

# 3 - TagCountToFastqPlugin

# Description:  Converts a master tagCount file containing all the tags of interest for your species/experiment (i.e., all of the tags with a minimum count greater than the -c parameter used in the MergeMultipleTagCountPlugin) from binary (.cnt) format into a FASTQ format file (.fq) that can then be used as input to one of the aligners BWA or bowtie2.
# Input: A binary tag count (.cnt) file containing all tags of interest (= master tag list).
# Output: The master tag list in FASTQ format (.fq). Can be used as input to BWA or bowtie2.

# 4 - SAMConverterPlugin

# Description:  Converts a SAM format alignment (.sam) file produced by one of the aligners, BWA or bowtie2, into a binary tagsOnPhysicalMap (.topm) file that can be used by the TagsToSNPByAlignmentPlugin for calling SNPs.
# Input: SAM format alignment (.sam) file produced by BWA or by the default (-M) mode of bowtie2.
# Output: A Binary tagsOnPhysicalMap (.topm) file that can be used by the TagsToSNPByAlignmentPlugin for calling SNPs.

# 5 - FastqToTBTPlugin:

# Description:  Generates a TagsByTaxa file for each FASTQ file in the input directory. One TagsByTaxa file is produced per FASTQ file. Requires a master list of tags of interest, which may come either from a tagCount (.cnt) or tagsOnPhysicalMap (.topm) file.
# Input: (1) FASTQ directory, (2) barcode key file, and (3) the master tag list in the form of either a binary tagCount (.cnt) file or a tagsOnPhysicalMap (.topm) file
# Output: Directory containing a corresponding tagsByTaxa file for every FASTQ file in the input directory

# 6 - MergeTagsByTaxaFilesPlugin

# Description:  Merges all .tbt.bin and/or (preferably) .tbt.byte files present in the input directory and all of its subdirectories. For the best genotyping results (proper calling of heterozygotes), we recommend using .tbt.byte files as input (produced by the FastqToTBTPlugin using the -y option). Note: to merge replicate samples, invoke "-x" - although it is not recommended as sample mix-ups are possible.
# Input: Directory containing multiple tagsByTaxa (.tbt.byte or .tbt.bin) files (produced by FastqToTBTPlugin). 
# Output: Merged tagsByTaxa file (it is best to send this to a separate directory from the input directory).

# 7 - tbt2vcfPlugin

# Description: Calls genotypes calls from TBT and TOPM files, one chromosome at a time
# Input: (1) tagsByTaxa (.tbt.byte) and (2) TagsOnPhysicalMap (.topm) files
# Output: VCF file for each chromosome present in the reference genome

#########################################################################################################################
