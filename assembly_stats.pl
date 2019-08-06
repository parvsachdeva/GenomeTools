#!/usr/bin/perl
use Cwd qw(cwd);
##############################################################################
#Program to calculate number of contigs, genome size, N50 and L50
##############################################################################
sub extractFileNameFromPath{
  my $path=$_[0];
  @aa=split '/',$path;
  return $aa[scalar(@aa)-1];
}
#######################################Writing data for making plots#######################################
sub plot{
  my $dir=cwd;
  open $fhR, ">$dir/tmp_plot_contigs.r";
  $output_plot=extractFileNameFromPath($ARGV[0]);
  $output_plot=~s/(\.fasta)//g;
  $output_plot=~s/(\.fa)//g;
  print "$dir\n";
  print $fhR "contig_lengths=read.table('tmp_plots.txt')\npdf('$dir/$output_plot.pdf')\nhist(log2(contig_lengths[,1]),breaks=40, col='pink', xlab='Log2 of contig lengths', main='Histogram of Log2 of Contig Lengths')";
  open $fh2,">$dir/tmp_plots.txt";
  my @array=@{$_[0]};
  foreach my $x (@array){
    print $fh2 "$x\n";
  }
  #print "***Creating a temporary file for generating plot***\n";
  system("Rscript $dir/tmp_plot_contigs.r");
  #print "***Removing the temporary file for plot***\n";
  system("rm $dir/tmp_plots.txt");
  system("rm $dir/tmp_plot_contigs.r");
  print "Plot created : $dir/$output_plot.pdf\n";
}
#######################################Removing Contigs Less than a Cut-Off#######################################
sub rem_contigs{
  my @lens=@{$_[0]};
  my $cutoff=$_[1];
  my $i=0;
  foreach my $x (@lens){
    if($x>$cutoff){
      $cor_sorted_contig_len[$i]=$x;
      #if($i<10){print "$x\n";}
      $i++;
    }
  }
  return(@cor_sorted_contig_len);
}
#######################################Computing Mean Contig Length#######################################
sub mean{
  my $mean=0;
  @array = @{$_[0]};
  foreach $x (@array){
    $mean+=$x;
  }
  $mean/=scalar @array;
  return int $mean;
}
##############################################################################
sub sum{
  my @array=@{$_[0]};
  my $sum=0;
  foreach my $x (@array){
    $sum+=$x;
  }
  return $sum;
}
#######################################Computing N50 and L50#######################################
sub n50_l50{
  my $n50=0, $i=0, $y=0, $j=0;
  $l50=0;
  my @array=@{$_[0]};
  my $genome_size=$_[1];
  for($i=scalar @array-1;$i>-1;$i--){
    $add+=$array[$i];
    #N10 & L10
    if($add>=$genome_size/10 && !$f10){
      $n10=$array[$i];
      $l10=scalar @array-$i;
      $f10=1;
    }
    #N20 & L20
    if($add>=$genome_size/5 && !$f20){
      $n20=$array[$i];
      $l20=scalar @array-$i;
      $f20=1;
    }
    #N30 & L30
    if($add>=$genome_size/3.333 && !$f30){
      $n30=$array[$i];
      $l30=scalar @array-$i;
      $f30=1;
    }
    #N40 & L40
    if($add>=$genome_size/2.5 && !$f40){
      $n40=$array[$i];
      $l40=scalar @array-$i;
      $f40=1;
    }
    #N50 & L50
    if($add>=$genome_size/2 && !$f50){
      $n50=$array[$i];
      $l50=scalar @array-$i;
      $f50=1;
    }
    #N60 & L60
    if($add>=$genome_size*(6/10) && !$f60){
      $n60=$array[$i];
      $l60=scalar @array-$i;
      $f60=1;
    }
    #N70 & L70
    if($add>=$genome_size*(7/10) && !$f70){
      $n70=$array[$i];
      $l70=scalar @array-$i;
      $f70=1;
    }
    #N80 & L80
    if($add>=$genome_size*(8/10) && !$f80){
      $n80=$array[$i];
      $l80=scalar @array-$i;
      $f80=1;
    }
    #N90 & L90
    if($add>=$genome_size*(9/10) && !$f90){
      $n90=$array[$i];
      $l90=scalar @array-$i;
      $f90=1;
      last;
    }

  }
  print "N10 : $n10\tL10 : $l10\n";
  print "N20 : $n20\tL20 : $l20\n";
  print "N30 : $n30\tL30 : $l30\n";
  print "N40 : $n40\tL40 : $l40\n";
  print "N50 : $n50\tL50 : $l50\n";
  print "N60 : $n60\tL60 : $l60\n";
  print "N70 : $n70\tL70 : $l70\n";
  print "N80 : $n80\tL80 : $l80\n";
  print "N90 : $n90\tL90 : $l90\n";
  #my @ret_val=($n50,$l50);
  #return @ret_val;
}
#######################################Opening the file#######################################
$file=$ARGV[0];
open $fh, "<$file" or die print "Not enough parameters given $!";
$seq_cutoff=$ARGV[1];
#######################################Variable declarations#######################################
$contig_count=-1;
$genome_size=0;
$n50=0;$l50=0;
$total_n_count=0;
$gc_content=0;
$at_content=0;
$atgc_content=0;
#######################################Parsing the file for info#######################################
while(!eof($fh)){
  $buf=<$fh>;
  chomp $buf;
  if($buf=~"^>"){
    $contig_count++;
    $n_count[$contig_count]=0;
  }else{
    $len[$contig_count]+=length($buf);
    $genome_size+=length($buf);
    $n_count[$contig_count]+=$buf=~tr/N|n//;
    $at_content+=$buf=~tr/A|T|a|t//;
    $gc_content+=$buf=~tr/G|C|g|c//;
  }
}

################################Sorting contig/scaffold length array################################
@len2=sort{ $a <=> $b }(@len);
@len3=rem_contigs(\@len2,$seq_cutoff);
$assembly_size=sum(\@len3);
$total_contigs=scalar @len3;
@ret=n50_l50(\@len3,$assembly_size);
$n50=$ret[0];$l50=$ret[1];
#$l50=$total_contigs-$l50+1;
#$n50=$ret_val[0];$l50=$ret_val[1];

$mean_contig_length=mean(\@len3);

$total_n_count=sum(\@n_count);
$at_content/=$genome_size/100;
$gc_content/=$genome_size/100;
$atgc_content=$at_content+$gc_content;
$other_bases=100-$atgc_content;
$average_n_per_contig=$total_n_count/$contig_count;
#######################################Printing all results#######################################
print "\nNumber of contigs : $total_contigs\n";
print "Assembly size : $assembly_size\n";
#print "N50 : $n50\n";
#print "L50 : $l50\n";
print "Longest contig length : $len2[-1]\n";
print "Mean contig length : $mean_contig_length\n";
print "Total number of N's : $total_n_count\n";
print "Average number of N's per contig : $average_n_per_contig\n";
print "IUPAC bases other than ATGC : $other_bases %\n";
print "GC content : $gc_content %\n";
#######################################Potting histogram of contig lengths#######################################
plot(\@len);
#######################################END#######################################
