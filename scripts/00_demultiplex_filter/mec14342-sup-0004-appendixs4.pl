#!/usr/bin/perl -w
#GJR 11/19/2012
#apply filters to scattered .vcf files created by "ParallelGenotype.pl", which splits up GATK genotyping jobs and submits to PBS queue
#filters include: 1) overall snp quality, 2) minimum overall depth, 3) minimum frequency of rare variant, 4) bi-allelic loci
 # and finally 5) a filter to remove all loci where read counts for apparent heterozygotes deviate from a binomial distribution with p=0.5 
#note: should save output files ending in .vcf so that GATK will recognize them during the next merging step
#Reproductive isolation and environmental adaptation shape the phylogeography of mountain pine beetle (Dendroctonus ponderosae)
#Eddy J. Dowle1,2*, Ryan R. Bracewell3, Michael E. Pfrender4, Karen E. Mock5, Barbara J. Bentz5,6, Gregory J. Ragland1,2*
#1 Department of Entomology, Kansas State University
#2 Department of Integrative Biology, University of Colorado, Denver
#3 Department of Integrative Biology, University of California, Berkeley 
#4 Department of Biological Sciences, University of Notre Dame 
#5 Department of Wildland Resources, Utah State University, 
#6 USDA Forest Service, Rocky Mountain Research Station, Logan UT
# *Correspondence:
#Eddy Dowle, Department of Integrative Biology, University of Colorado Denver, 1151 Arapahoe, SI 2071 Denver, CO 80204, USA. email: edwina.dowle@ucdenver.edu
#Gregory Ragland, Department of Integrative Biology, University of Colorado, Denver, 1151 Arapahoe, SI 2071, Denver, CO 80204, USA. email: gregory.ragland@ucdenver.edu

# specify packages and specify filter variables
use strict;
#use lib "/afs/crc.nd.edu/user/g/gragland/PerlLibs/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
#perl package with cummulative frequency distributions (including binomial), install from CRAN
use Math::CDF;
#my %GoodContigs;
my $minAltFreq=0.05;
my $minDepth=50;
my $binomThresh=0.05;
my $minQual=21;


die "usage: perl FileterGatkVcfs.pl infile " unless @ARGV ==1;

my $infile=shift @ARGV;

#my $infile = "GATKgroup1.vcf";
#open INFILE, "<$infile" or die $1;


#read directory and store .vcf files to be filtered
#my @dir;
#opendir (DIR, ".") or die $!;
#while (my $file = readdir DIR) {push @dir, $file};
#closedir(DIR);
#@dir = grep /.vcf\.\d+$/, @dir;

for ($infile) {

   open INFILE, "<$infile";
   my $outfile = $infile;
   $outfile =~ s/\.vcf//; ##???
#  open OUTFILE, my $outfile = ">$infile.filtered.vcf";
   open OUTFILE, ">$outfile.filtered.vcf";
  while (<INFILE>) {
    if (/^\#/) {
      print OUTFILE "$_";
      next;
    }
    chomp;
    my @vals=split "\t";

    #primary filters for minimum depth, minimum quality, bi-allelic loci, and minimum alternate allele frequency
    #filters entries with multiple alternate alleles
    next if $vals[4] =~ m/\S\S/;
    my $qual=$vals[5];
    next if $qual < $minQual;
    /DP=(\d+)/;
    next if $1 < $minDepth;
    /AF=([\d\.e\-]+)/;
    next if $1 < $minAltFreq and $1 < (1-$minAltFreq);

    my $contig=$vals[0];
    my $allele1=0;
    my $allele2=0;
    for my $info (@vals[9..$#vals]) {
      next if $info =~ m/\.\/\./;
      my @AllInfo=split ":", $info;
      my @Counts=split ",", $AllInfo[1];
      if ($AllInfo[1] !~ m/0/) {
	$allele1=$allele1+$Counts[0];
	$allele2=$allele2+$Counts[1];
      }
    }
    my @acounts=sort {$a <=> $b} ($allele1,$allele2);
    my $prob = 0.5;
    if ($allele1 > 0 and $allele2 > 0) {
      $prob = &Math::CDF::pbinom($acounts[0], ($acounts[0]+$acounts[1]), 0.5);
    }
    if ( $prob < (1-$binomThresh) and $prob > $binomThresh) {
      print OUTFILE "$_\n";
    }
  }
  close INFILE;
  close OUTFILE;

}