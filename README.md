# GenomeTools
Tools to analyse and extract genomic data.

## 1. assembly_stats.pl :
Perl and R script to get assembly/ fasta info. The R script is embedded in the Perl code and is used for plotting the histogram of contig lengths.
##### Syntax : ./assembly_stats.pl your_assembly.fasta min_contig_length_to_consider
##### Output : 
L10-L90, N10-N90, histogram of contig/ scaffold lengths(need an installation of R), number of contigs, assembly size, longest contig length, mean contig length, total number of N's, average number of N's per contig, IUPAC bases other than ATGC(%), GC content(%).


## 2. get_scaffold.py
Python script to fetch a specific chromosome/ scaffold/ contig from a fasta file. Use Python v3 or above.
##### Syntax : python get_scaffold.py -i input.fasta -s scaffold_or_chromosome_name -o output.fasta
##### Output : 
Fasta file with specified chromosome/ scaffold header and sequence.

## 3. rev_com.py
Python script to get complement, reverse or reverse complement of sequences in a fasta file. The masked sequences will not be affected and the case is maintained. Use Python v3 or above.
##### Syntax : python rev_com.py -i input.fasta -c choice -o output.fasta
(-c : 0 for reverse; 1 for complement; 2 for reverse complement)
##### Output : 
Fasta file with desired operation performed on all the sequences in the file.

##  4. variant_calling_pipeline.sh
Shell script for complete variant calling process. This script requires path to GATK (jar file), PICARD (jar file) and a global installation of BWA (http://bio-bwa.sourceforge.net/). It will ask for an input of the names of the paired-end read files (fastq) and the genome file (fasta).
The following steps will be carried out : 
  1. Genome indexing
  2. Quality control and trimming of reads
  3. Alignment of reads against the genome
  4. Sorting SAM file by coordinate and conversion to BAM
  5. Getting sequence depth
  6. Building index for BAM file
  7. Creating realignment targets
  8. Realigning indels
  9. Variant calling
  10. Extracting SNPs and INDELS
  11. Filtering SNPs and INDELS
  12. Base Quality Score Recalibration (BQSR)
  13. Calling final variants
  14. Extracting final SNPs and INDELS
  15. Final filtering of SNPs and INDELS
 
 #### Syntax : bash variant_calling_pipeline.sh
 #### Output : 
VCF files containing SNPs and INDELS present in sample w.r.t. the genome.

## 5. syn_gen.pl
Perl script to create a synthetic genome/ chromosome using a random number generator.
The program will ask for the number of bases required and a file name to write the output to. Only ATGC bases are printed in the file in random.
#### Syntax : ./syn_gen.py
#### Output : 
A fasta file with a header consisting of a randomly generated genome.

## 6. shortest_common_superstring.pl
Perl script that can :
  1. Find all possible k-mers in a given sequence
  2. Find all unique k-mers
  3. Find the shortest common superstring using the unique k-mers
 The program will ask for a string and the size of k-mers to compute the above.
 #### Syntax : ./shortest_common_superstring.pl
 #### Output : 
 Unique k-mers, Shortest Common Superstring
 
