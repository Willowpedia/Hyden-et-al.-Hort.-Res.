#########################################################################################################################
#########################################################################################################################
#
# Project: Willow GBS Pipeline script
# Reference genome: Female 'Salix purpurea' clone 94006 v4 (PacBio polished assembly)
# Population: 'Association and F2 family'
# Date: 2018-03-18
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

# Modified by Brennan Hyden 11.5.2018
#
#########################################################################################################################

#########################################################################################################################
#################################################### Begin GBS Pipeline #################################################

# Set your working directory
cd "/Volumes/HAL9000/genomics/salix_ngs_fastq/fastq_gbs"

####################################################
# Make directories and index reference genome
####################################################

# Remember to use ">chr1, >chr2 ... >chr21" FASTA header formats, TASSEL will not accept any modifications
# Rename chromosomes so they're in order (chr1 -> chr...N, where N=number of sequences)
awk '/^>/{print ">chr" ++i; next}{print}' < genome/Salix_purpurea_var_94006_mainGenome.fasta > genome/Salix_purpurea_var_94006_v4.mainGenome.fasta

# Print fasta ids of the original and modified referece
grep -o -E "^>\w+" genome/Salix_purpurea_var_94006_mainGenome.fasta | tr -d ">" > genome/Salix_94006v4_mainGenome.chrIDs_original.txt
grep -o -E "^>\w+" genome/Salix_purpurea_var_94006_v4.mainGenome.fasta | tr -d ">" > genome/Salix_purpurea_var_94006_v4.mainGenome.chrIDs_modified.txt

# Concatenate reference seq id files for later use
paste genome/Salix_94006v4_mainGenome.chrIDs_original.txt genome/Salix_purpurea_var_94006_v4.mainGenome.chrIDs_modified.txt > genome/Salix_94006_v4.mainGenome.chrIDs_orig_mod.txt

####################################################
# Index reference genome
####################################################

# Normalize fasta
picard NormalizeFasta I=Reference_genome/Fv5_main_genome_chr15W.fasta O=Reference_genome/Fv5_mainGenome_chr15W.fasta

# Build a BWA-index for the reference genome
bwa index -a is Reference_genome/Fv5_mainGenome_chr15W.fasta

# Build a SAMTools index for the reference genome
samtools faidx Reference_genome/Fv5_mainGenome_chr15W.fasta

####################################################
# Get mergedTagCount files from raw fastq
####################################################

cd "Assoc_F2_ApeKI"

# FastqToTagCountPlugin 
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -FastqToTagCountPlugin -i fastq -k AoB.keyfile -e ApeKI -s 350000000 -c 1 -o tagCounts -endPlugin -runfork1 > tagCounts/FastqToTagCountPlugin.log 2> tagCounts/FastqToTagCountPlugin.err

# MergeMultipleTagCountPlugin 
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeMultipleTagCountPlugin -i tagCounts -o mergedTagCounts/mergedTagCounts.cnt -c 5 -endPlugin -runfork1 > mergedTagCounts/mergedTagCounts_c5_nb.log 2> mergedTagCounts/mergedTagCounts_c5_nb.err

# BinaryToTextPlugin
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -BinaryToTextPlugin -i mergedTagCounts/mergedTagCounts.cnt -o mergedTagCounts/mergedTagCounts.cnt.txt -t tagCounts -endPlugin -runfork1

# FastqToTBTPlugin
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -FastqToTBTPlugin -i fastq -y -k AoB.keyfile -e ApeKI -o tbt -t mergedTagCounts/mergedTagCounts.cnt -endPlugin -runfork1 > tbt/tbt0.log 2> tbt/tbt0.err

# MergeTagsByTaxaFilesPlugin
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeTagsByTaxaFilesPlugin -i tbt -s 350000000 -o mergedTBT/mergedTBT.tbt.byte -endPlugin -runfork1 > mergedTBT/tbt_nb.log 2> mergedTBT/tbt_nb.err

# MergeTagsByTaxaFilesPlugin (-x merge tag counts of taxa with identical short names)
# ~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -MergeTagsByTaxaFilesPlugin -i tbt -x -s 350000000 -o mergedTBT/mergedTBT.x.tbt.byte -endPlugin -runfork1 > mergedTBT/tbt_x_nb.log 2> mergedTBT/tbt_x_nb.err

# TagCountToFastqPlugin
~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx56g -fork1 -TagCountToFastqPlugin -i mergedTagCounts/mergedTagCounts.cnt -o mergedTagCounts/mergedTagCounts.fq -c 5 -endPlugin -runfork1 > mergedTagCounts/mergedTagCounts_c5_cvtnb.log 2> mergedTagCounts/mergedTagCounts_c5_cvtnb.err

# BinaryToTextPlugin
# ~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx64g -fork1 -BinaryToTextPlugin -i mergedTBT/mergedTBT.x.tbt.byte -o mergedTBT/AoB_94006v4.mergedTBT.txt -t TBTByte -endPlugin -runfork1
# gzip -c mergedTBT/AoB_94006v4.mergedTBT.txt > mergedTBT/AoB_94006v4.mergedTBT.txt.gz

# ~/programs/tassel3-standalone/run_pipeline.pl -Xms2g -Xmx64g -fork1 -BinaryToTextPlugin -i mergedTBT/mergedTBT.x.tbt.byte -o mergedTBT/AoB_94006v4.EcoT22I.mergedTBT.txt -t TBTByte -endPlugin -runfork1
# gzip -c mergedTBT/AoB_94006v4.EcoT22I.mergedTBT.txt > mergedTBT/AoB_94006v4.EcoT22I.mergedTBT.txt.gz

####################################################
# Map mergedTagCounts to reference (-o SAM)
####################################################

cd "Hyden/Assoc_F2_ApeKI"

####################
# BWA MEM
####################

# BWA-mem Perl Script (-i Fastq [.fq] -o SequenceAlignment [.sam])

bwa mem -t 12 Reference_genome/Fv5_mainGenome_chr15W.fasta mergedTagCounts/mergedTagCounts.fq > topm/Fv5_mainGenome_chr15W.sam

# First run CHMOD on that Perl shiznit ("chmod +x" allows us to execute the script)!
chmod +x filter_bwa_mem.pl

# Convert from BWA-format to a readable Bowtie2-format SAM file
perl ./filter_bwa_mem.pl topm/Fv5_mainGenome_chr15W.sam topm/Fv5_mainGenome_chr15W.filt.sam 1

# SAMConverterPlugin - BWA
~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -SAMConverterPlugin -i topm/Fv5_mainGenome_chr15W.filt.sam -o topm/Fv5_mainGenome_chr15W.topm -endPlugin -runfork1 > topm/topm_c5.log 2> topm/topm_c5.err

####################################################
# "Taxa UNMERGED" - Call variants (-o VCF) 
####################################################

# tbt2vcfPlugin: "UNMERGED" TBT (-i mergedTBT/mergedTBT.tbt.byte)
~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -tbt2vcfPlugin -i mergedTBT/mergedTBT.tbt.byte -m topm/Fv5_mainGenome_chr15W.topm -o tbt2vcf -mnMAF 0.05 -mnLCov 0.0 -ak 4 -s 1 -e 844 -endPlugin -runfork1 > tbt2vcf/tbt2vcf.log 2> tbt2vcf/tbt2vcf.err

# Set working directory to tbt2vcf output
cd "Hyden/Assoc_F2_ApeKI"

# Compress the mergedTBT for each chromosome
gzip tbt2vcf/mergedTBT.c*

# Give vcftools the path to PERL5
export PERL5LIB=~/programs/vcftools_0.1.13/lib/perl5/site_perl

# Concatenate VCFs and gzip the output
~/programs/vcftools_0.1.13/bin/vcf-concat tbt2vcf/mergedTBT.c*.gz | gzip -c > hyden_vcf/Fv5_mainGenome_chr15W.ape.vcf.gz

# De-compress VCF to merge and sort
gunzip -cd hyden_vcf/Fv5_mainGenome_chr15W.ape.vcf.gz > hyden_vcf/Fv5_mainGenome_chr15W.ape.vcf

# Merge duplicate SNPs
~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -MergeDuplicateSNP_vcf_Plugin -i hyden_vcf/Fv5_mainGenome_chr15W.ape.vcf -ak 4 -o hyden_vcf/Fv5_mainGenome_chr15W.ape.MergedDup.vcf -endPlugin -runfork1

# Sort the VCF
~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -SortGenotypeFilePlugin -inputFile hyden_vcf/Fv5_mainGenome_chr15W.ape.MergedDup.vcf -outputFile hyden_vcf/Fv5_mainGenome_chr15W.ape.MergedDupSorted.vcf -fileType VCF -endPlugin -runfork1


# Compress and index sorted VCF
bgzip hyden_vcf/Fv5_mainGenome_chr15W.ape.MergedDupSorted.vcf
tabix -p vcf hyden_vcf/Fv5_mainGenome_chr15W.ape.MergedDupSorted.vcf.gz

# TASSEL GUI

~/programs/tassel-5-standalone/start_tassel.pl -Xms4g -Xmx56g

# Family 317 (maxMAF=0.85, minMAF=0.15, miss=200, maxHet=0.85, minHet=0.15)
#bgzip vcf/Fv5_mainGenome_chr15W.family317.filt.vcf
#tabix -p vcf vcf/Fv5_mainGenome_chr15W.mG.family317.filt.vcf.gz

# Association panel (maxMAF=1, minMAF=0.05, miss=80)
#bgzip vcf/Fv5_mainGenome_chr15W.mG.assoc.filt.vcf
#tabix -p vcf vcf/Fv5_mainGenome_chr15W.mG.assoc.filt.vcf.gz

####################################################
# "Taxa MERGED" - Call variants (-o VCF) 
####################################################

# tbt2vcfPlugin: "MERGED" TBTx (-i mergedTBT/mergedTBTx.tbt.byte)
#~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -tbt2vcfPlugin -i mergedTBT/mergedTBT.tbt.byte -m topm/Fv5_mainGenome_chr15W.topm -o tbtx2vcf -mnMAF 0.05 -mnLCov 0.0 -ak 4 -s 1 -e 844 -endPlugin -runfork1 > tbtx2vcf/tbtx2vcf.log 2> tbtx2vcf/tbtx2vcf.err

# Set working directory to tbt2vcf output
#cd "Hyden/Assoc_F2_ApeKI"

# Compress the mergedTBT for each chromosome
#gzip tbtx2vcf/mergedTBT.c*

# Give vcftools the path to PERL v5
#export PERL5LIB=~/programs/vcftools_0.1.13/lib/perl5/site_perl

# Concatenate VCFs and gzip the output
#~/programs/vcftools_0.1.13/bin/vcf-concat tbtx2vcf/mergedTBT.c*.gz | gzip -c > vcf_mergedTaxa/Fv5_mainGenome_chr15W.vcf.gz

# De-compress VCF to merge and sort
#gunzip -cd vcf_mergedTaxa/Fv5_mainGenome_chr15W.vcf.gz > vcf_mergedTaxa/Fv5_mainGenome_chr15W.vcf

# Merge duplicate SNPs
#~/programs/tassel3-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -MergeDuplicateSNP_vcf_Plugin -i vcf_mergedTaxa/Fv5_mainGenome_chr15W.vcf -ak 4 -o vcf_mergedTaxa/Fv5_mainGenome_chr15W.mergedDup.vcf -endPlugin -runfork1

# Sort the VCF
#~/programs/tassel-5-standalone/run_pipeline.pl -Xms4g -Xmx56g -fork1 -SortGenotypeFilePlugin -inputFile vcf_mergedTaxa/Fv5_mainGenome_chr15W.mergedDup.vcf -outputFile vcf_mergedTaxa/AFv5_mainGenome_chr15W.mergedDupSorted.vcf -fileType VCF -endPlugin -runfork1

# Compress and index sorted VCF
#bgzip vcf_mergedTaxa/Fv5_mainGenome_chr15W.mergedDupSorted.vcf
#tabix -p vcf vcf_mergedTaxa/Fv5_mainGenome_chr15W.mergedDupSorted.vcf.gz

# TASSEL GUI
#~/programs/tassel-5-standalone/start_tassel.pl -Xms4g -Xmx56g

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

#bgzip Spurpurea94006v5.1.gene_exons.gff3
#tabix -s 1 -b 2 -e 3 Spurpurea94006v5.1.gene_exons.gff.gz

#bcftools annotate -a Spurpurea94006v5.1.gene_exons.gff.gz -c CHR,START,END,GENE -O z -h key=INFO,ID=ANN,Number=1,Type=Integer,Description='Salix purpurea var 94006 v5 gene annotation' --threads 12 -o Fv5_main_genome_chr15W_mergedDupSorted.gff.vcf Fv5_main_genome_chr15W_mergedDupSorted.vcf

####################################################
# Accessory Hapmap functions                     
####################################################

# TagsToSNPByAlignmentPlugin (-i [.tbt.byte] -o [.hmp] -m [.topm] -mUpd [withVariants.topm])
# ~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -TagsToSNPByAlignmentPlugin -i mergedTBT/mergedTBT.tbt.byte -ref genome/Salix_94006v4_mainGenome.fasta -m topm/94006_v2.topm -y -mnMAF 0.0 -inclGaps -callBiSNPsWGap -sC 1 -eC 246 -o hapmap/Spur_94006_v3_chr+.hmp.txt -endPlugin -runfork1

# MergeIdenticalTaxaPlugin
#~/programs/tassel3-standalone/run_pipeline.pl -Xmx64g -fork1 -MergeIdenticalTaxaPlugin -hmp vcf/AoB_94006_v4.mG.assoc.filt.hmp.txt -o AoB_94006_v4.mG.assoc.mergedTaxa.hmp.txt -hetFreq 0.75 -sC 1 -eC 363 -endPlugin -runfork1

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
