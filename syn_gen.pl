#!/usr/bin/perl

print "Enter size of genome : ";
$gsize=<>;
print "Enter file to write output : ";
$fname=<>;
$gsize=10000;
$fname="syn_gen.fa";
@bases=("A","T","G","C");
for(my $i=0;$i<$gsize;$i++){
  my $rb=int rand(4);
  $genome.=$bases[$rb];
}

open($fh1,">$fname");
print $fh1 ">Synthetic_Genome_length_$gsize\n";
print $fh1 $genome;
