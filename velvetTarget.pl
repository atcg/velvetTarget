#!/usr/bin/perl

# The MIT License (MIT)
# 
# Copyright (c) 2013 Evan Melstad (evanmelstad@ucla.edu)
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $miSeq = 0;
my $name;
my $R1file;
my $R2file;
my $singlesFile;
my $fromK = 19; #default starting Kmer value is 19
my $toK = 201; #default ending emer value is 201
my $phredBase = 33; #defaults to phred33 quality scores
my $outFile;
my $usage = "Usage: perl velvetTarget.pl -i <infile.txt> -o <outfile.txt>\n";

GetOptions  ("name=s"    => \$name,
             "R1=s"      => \$R1file,
             "R2=s"      => \$R2file,
             "singles=s" => \$singlesFile,
             "from=i"    => \$fromK,
             "to=i"      => \$toK,
             "phred=i"   => \$phredBase,
             "out=s"     => \$outFile,
             "miseq"     => \$miSeq,
             "help|man"  => \$help,) || pod2usage(2);

pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if $help;

if ((!$name) or (!$R1file) or (!$R2file) or (!$singlesFile) or (!$fromK) or (!$toK) or (!$outFile)) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}
if (($fromK % 2 != 1) or ($toK % 2 != 1)) {
    print "\n***Kmer value bounds (--from and --to) must both be odd integers.***\n\n";
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}




# Run velveth for all odd numbers between the specified kmer values
my @kmer_choices = grep {$_ % 2 == 1} $fromK..$toK;

foreach my $kmer (@kmer_choices) {
    my $velvet_dir_name = $name . "_" . $kmer . "_velvet";
    my $dbName = $name . "_" . $kmer;
    if ($miSeq) {
        system("velveth $velvet_dir_name $kmer -short -fastq $singlesFile -longPaired -separate -fastq $R1file $R2file");
        system("velvetg $velvet_dir_name -exp_cov auto -cov_cutoff auto");
    } else {
        system("velveth $velvet_dir_name $kmer -short -fastq $singlesFile -shortPaired -separate -fastq $R1file $R2file");
        system("velvetg $velvet_dir_name -exp_cov auto -cov_cutoff auto");
    }
    system("makeblastdb -in $velvet_dir_name/contigs.fa -dbtype nucl -out $velvet_dir_name/$dbName");        
}






























#Documentation
__END__

=head1 NAME
velvetTarget.pl

=head1 SYNOPSIS 

perl velvetTarget.pl --R1 <file> --R2 <file> --singles <file> --from \
<int> --to <int> --phred <int> --out <file>

 Options:
   -R1=s            R1 reads file
   -R2=s            R2 reads file
   -singles=s       Singletons and joined reads file
   -from=i          Starting kmer value (default 19)
   -to=i            Ending kmer value (default 201)
   -phred=i         Phred quality scoring (default = 33)
   -out=s           Log file name
   -miseq           Include this flag if you have MiSeq reads over 200bp
   -help|man        documentation


=head1 DESCRIPTION

Choosing the right kmer value is critical for de novo assembly of
contigs from next generation sequencing, but this is difficult. This
task is easier for genome reconstruction, as the goal is usually to
assemble the longest contigs possible (VelvetOptimiser.pl is great for
this). For target enrichment experiments, however, the goal is rather to
recover as many of the targets as completely as possible.

This script attempts to find the k value that achieves this. Instead of
basing the choice of k off of summary statistics of the assembled
contigs themselves, it rather blasts the baits or target regions against
the assembly, and returns the k value and assembly that recovers the
most targets the most completely. 

=cut