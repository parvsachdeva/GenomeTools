#!/bin/bash
##PreRequisites
#wd should have 2 .fastq.gz files (raw reads) and 1 .fa file (reference)
#****Enter complete path to picard.jar****
picard="/home/ibab/Bioinformatics_Softwares/picard/picard.jar"
#****Enter path to GATK.jar****
gatk="/home/ibab/Bioinformatics_Softwares/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"


##Asking User for File Names and Moving to Respective Directories
echo "****Enter first raw read data file name (named as *.fastq)****"
read raw_file_1
echo "****Enter second raw read data file name (named as *.fastq)****"
read raw_file_2
echo "****Enter reference genome file name (named as *.fa)****"
read ref_file
echo "****Enter Type of Raw Data (e.g. cancer/disease)****"
read dtype

##Making Directories to Order Work
mkdir rawdata
mkdir reference
mkdir alignment
mkdir trimmed
mkdir variants
raw=$PWD/rawdata
ref=$PWD/reference
trim=$PWD/trimmed
align=$PWD/alignment
var=$PWD/variants
mv $raw_file_1 $raw
mv $raw_file_2 $raw
mv $ref_file $ref

##Indexing the Genome
echo "****BWA INDEXING REFERENCE GENOME****"
echo "*****-2*****"
bwa index $ref/$ref_file
echo "*****-1*****"
samtools faidx $ref/$ref_file
ref_file_dict=${ref_file/".fa"/".dict"}
echo "*****0*****"
java -jar $picard CreateSequenceDictionary R=$ref/$ref_file O=$ref/$ref_file_dict

##Quality Control and Trimming
fastqc $raw/*.gz
trim_galore --fastqc --paired $raw/$raw_file_1 $raw/$raw_file_2 -o $trim/

trim_file_1=${raw_file_1/".fastq"/"_val_1.fq"}
trim_file_2=${raw_file_2/".fastq"/"_val_2.fq"}

##Aligning Against the Genome
echo "*****1*****"
bwa mem -M -R "@RG\\tID:cancer\\tLB:cancer\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:cancer" $ref/$ref_file $trim/$trim_file_1 $trim/$trim_file_2 > $align/$dtype".sam"

##Sorting SAM File by Coordinate and converting to BAM
echo "*****2*****"
java -jar $picard SortSam INPUT=$align/$dtype".sam" OUTPUT=$align/$dtype"_sorted.bam" SORT_ORDER=coordinate

##Getting Sequence Depth
echo "*****3*****"
samtools depth -a $align/$dtype"_sorted.bam" > $align/$dtype"_depth.txt"
echo "*****4*****"
java -jar $picard MarkDuplicates INPUT=$align/$dtype"_sorted.bam" OUTPUT= $align/$dtype"_sorted_dedup.bam" METRICS_FILE=$align/$dtype"_metrics.txt"

##Building Index for BAM File
echo "*****5*****"
java -jar $picard BuildBamIndex INPUT= $align/$dtype"_sorted_dedup.bam"

##Creating Realignment Targets
echo "*****6*****"
java -jar $gatk \
-T RealignerTargetCreator \
-R $ref/$ref_file \
-I $align/$dtype"_sorted_dedup.bam" \
-o $align/$dtype"_realignment_targets.list"

##Realigning Indels
echo "*****7*****"
java -jar $gatk \
-T IndelRealigner \
-R $ref/$ref_file \
-I $align/$dtype"_sorted_dedup.bam" \
-targetIntervals $align/$dtype"_realignment_targets.list" \
-o $align/$dtype"_realigned_reads.bam"

##Variant Calling
echo "*****8*****"
java -jar $gatk \
-T HaplotypeCaller \
-R $ref/$ref_file \
-I $align/$dtype"_realigned_reads.bam" \
-o $var/$dtype"_raw_variants.vcf"

##Extracting SNPs and INDELS
#extracting SNPs
echo "*****9*****"
java -jar $gatk \
-T SelectVariants \
-R $ref/$ref_file \
-V $var/$dtype"_raw_variants.vcf" \
-selectType SNP \
-o $var/$dtype"_raw_snps.vcf"

#extracting INDELs
echo "*****10*****"
java -jar $gatk \
-T SelectVariants \
-R $ref/$ref_file \
-V $var/$dtype"_raw_variants.vcf" \
-selectType INDEL \
-o $var/$dtype"_raw_indels.vcf"

##Filtering SNPs and INDELs
#filtering SNPs
echo "*****11*****"
java -jar $gatk \
-T VariantFiltration \
-R $ref/$ref_file \
-V $var/$dtype"_raw_snps.vcf" \
--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
--filterName "basic_snp_filter" \
-o $var/$dtype"_filtered_snps.vcf"

#filtering INDELs
echo "*****12*****"
java -jar $gatk \
-T VariantFiltration \
-R $ref/$ref_file \
-V $var/$dtype"_raw_indels.vcf" \
--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
--filterName "basic_indel_filter" \
-o $var/$dtype"_filtered_indels.vcf"

##Base Quality Score Recalibration (BQSR)
echo "*****13*****"
java -jar /$gatk \
-T BaseRecalibrator \
-R $ref/$ref_file \
-I $align/$dtype"_realigned_reads.bam" \
-knownSites $var/$dtype"_filtered_snps.vcf" \
-knownSites $var/$dtype"_filtered_indels.vcf" \
-o $var/$dtype"_recal_data.table"
echo "*****14*****"
java -jar $gatk \
-T PrintReads \
-R $ref/$ref_file \
-I $align/$dtype"_realigned_reads.bam" \
-BQSR $var/$dtype"_recal_data.table" \
-o $var/$dtype"_recal_reads.bam"
echo "*****15*****"
##Calling Final Variants
java -jar $gatk \
-T HaplotypeCaller \
-R $ref/$ref_file \
-I $var/$dtype"_recal_reads.bam" \
-o $var/$dtype"_raw_variants_recal.vcf"

##Extracting Final SNPs and INDELs
#extracting SNPs
echo "*****16*****"
java -jar $gatk \
-T SelectVariants \
-R $ref/$ref_file \
-V $var/$dtype"_raw_variants_recal.vcf" \
-selectType SNP \
-o $var/$dtype"_raw_snps_recal.vcf"

#extracting INDELs
echo "*****17*****"
java -jar $gatk \
-T SelectVariants \
-R $ref/$ref_file \
-V $var/$dtype"_raw_variants_recal.vcf" \
-selectType INDEL \
-o $var/$dtype"_raw_indels_recal.vcf"

##Final Filtering
#filtering SNPs
echo "*****18*****"
java -jar $gatk \
-T VariantFiltration \
-R $ref/$ref_file \
-V $var/$dtype"_raw_snps_recal.vcf" \
--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
--filterName "basic_snp_filter" \
-o $var/$dtype"_filtered_snps_final.vcf"

#filtering INDELs
echo "*****19*****"
java -jar $gatk \
-T VariantFiltration \
-R $ref/$ref_file \
-V $var/$dtype"_raw_indels_recal.vcf" \
--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
--filterName "basic_indel_filter" \
-o $var/$dtype"_filtered_indels_final.vcf"
