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
