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
use File::Basename qw(basename);
use Bio::SearchIO;
use Statistics::R;

my $help = 0;
my $miSeq = 0;
my $name;
my $adapterFile;
my $R1file;
my $R2file;
my $fromK = 19; #default starting Kmer value is 19
my $toK = 201; #default ending emer value is 201
my $probes;

#Print out command used to run script:
print "perl " . $0 . " ";foreach my $Argy (@ARGV) {print $Argy . " ";}print "\n\n\n";

GetOptions  ("name=s"           => \$name,
             "adapters=s"       => \$adapterFile,
             "R1=s"             => \$R1file,
             "R2=s"             => \$R2file,
             "from=i"           => \$fromK,
             "to=i"             => \$toK,
             "probes|targets=s" => \$probes,
             "miseq"            => \$miSeq,
             "help|man"         => \$help,) || pod2usage(2);

pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if $help;

if ((!$name) or (!$adapterFile) or (!$R1file) or (!$R2file) or (!$fromK)
    or (!$toK) or (!$probes)) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}
if (($fromK % 2 != 1) or ($toK % 2 != 1)) {
    print "\n***Kmer value bounds (--from and --to) must both be odd integers.***\n\n";
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my $scytheR1out = $name . ".R1.scythe";
my $scytheR2out = $name . ".R2.scythe";
system("scythe -a $adapterFile $R1file -q sanger -o $scytheR1out");
system("scythe -a $adapterFile $R2file -q sanger -o $scytheR2out");

my $sickleR1out = $name . ".R1.sickle";
my $sickleR2out = $name . ".R2.sickle";
my $sickleIndOut = $name . ".singles.sickle";
system("sickle pe -f $scytheR1out -r $scytheR2out -t sanger -o $sickleR1out -p $sickleR2out -s $sickleIndOut");
unlink $scytheR1out;
unlink $scytheR2out;

# Fastq-join
my $joined = $name . ".fqj.join.fastq";
my $R1_postjoin = $name . ".fqj.un1.fastq";
my $R2_postjoin = $name . ".fqj.un2.fastq";
system("fastq-join -v ' ' -m 10 $sickleR1out $sickleR2out -o $name.fqj.%.fastq");

my $singlesFile = $name . ".singles_and_joined.clean.fastq"; #created by concatenating
                                                             #sickle singletons and
                                                             #merged reads from fastq-join
system("cat $joined $sickleIndOut > $singlesFile");
unlink $sickleR1out;
unlink $sickleR2out;
unlink $sickleIndOut;
unlink $joined;

# Run velveth for all odd numbers between the specified kmer values
my @kmer_choices = grep {$_ % 2 == 1} $fromK..$toK;
my %blastingStats; # Holds information regarding the number of targets matching
                   # contigs and completeness of match
my @targets_with_hits_for_R;
my @targets_with_one_hit_for_R;
my @targets_with_one_HSP_for_R;
my @targets_nested_in_contig_for_R;

foreach my $kmer (@kmer_choices) {
    my $velvet_dir_name = $name . "_" . $kmer . "_velvet";
    my $dbName = $name . "_" . $kmer;
    my ($probeFilename) = basename($probes);
    my $blastOutName = $probeFilename . "_blasted_to_" . $dbName;
    my $blastStatsLog = $velvet_dir_name . "/$blastOutName" . "_blastStats.txt";
    
    if ($miSeq) {
        system("velveth $velvet_dir_name $kmer -short -fastq $singlesFile -longPaired -separate -fastq $R1_postjoin $R2_postjoin");
        system("velvetg $velvet_dir_name -exp_cov auto -cov_cutoff auto");
    } else {
        system("velveth $velvet_dir_name $kmer -short -fastq $singlesFile -shortPaired -separate -fastq $R1_postjoin $R2_postjoin");
        system("velvetg $velvet_dir_name -exp_cov auto -cov_cutoff auto");
    }

    system("makeblastdb -in $velvet_dir_name/contigs.fa -dbtype nucl -out $velvet_dir_name/$dbName");
    
    # Blast targets against assemblies
    system("blastn -db $velvet_dir_name/$dbName -query $probes -out $velvet_dir_name/$blastOutName");
    
    #Compile summary stats of blast results. We want to know:
    # 1. How many targets have at least one matching contig when blasted against the assembly?
    $blastingStats{'targets_with_hits'} = 0;
    # 2. How many targets have exactly one matching contig (hit) when blasted against the assembly?
    $blastingStats{'targets_with_one_hit'} = 0;
    # 3. How many targets with a single hit have exactly one HSP when blasted against the assembly?
    $blastingStats{'targets_with_one_hsp'} = 0;
    # 4. How many targets with a single hsp are matched to at least 98% of their
    # length (i.e. are nested within a contig with overlap)? Note this measure
    # might need tweaking. For instance, if the first basepair of the query sequence
    # is mismatched, but all remaining basepairs are a match, does blast chop off
    # the first basepair and report the alignment as, say 119bp out of a 120bp query.
    # This could be misleading if the first or last base pairs are simply mutations
    # from the query.
    $blastingStats{'targets_nested_within_contig'} = 0;

    my $blastIO = Bio::SearchIO->new(-format => 'blast',
                                     -file => "$velvet_dir_name/$blastOutName");
    while (my $blastResult = $blastIO->next_result()) {
        if ($blastResult->num_hits() > 0) {
            $blastingStats{'targets_with_hits'}++;
            if ($blastResult->num_hits() == 1) {
                $blastingStats{'targets_with_one_hit'}++;
                while (my $hit = $blastResult->next_hit()) {
                    if ($hit->num_hsps() == 1) {
                        $blastingStats{'targets_with_one_hsp'}++;
                        my $hsp = $hit->next_hsp();
                        if ((($hsp->length('total') / $blastResult->query_length()) > 0.98)
                            and ($hsp->start('subject') != 1)
                            and ($hsp->end('subject') != 1)) { #try changing
                                                               #'total' to 'hit'
                                                               #in the length call
                                                               #to see what happens
                            $blastingStats{'targets_nested_within_contig'}++
                        }
                    }            
                }       
            }   
        }                  
    }
    open(my $STATSFH, ">", $blastStatsLog) or die "Couldn't create the file to store the blast stats.\n";
    print $STATSFH "Number of target regions: " . $blastIO->result_count() . ".\n";
    print $STATSFH "Number of target regions with hits: $blastingStats{'targets_with_hits'}. ("
        . ((($blastingStats{'targets_with_hits'})/($blastIO->result_count))*100) . "% of total)\n";
    print $STATSFH "Number of target regions with ONE hit: $blastingStats{'targets_with_one_hit'} ("
        . ((($blastingStats{'targets_with_one_hit'})/($blastIO->result_count))*100) . "% of total)\n";
    print $STATSFH "Number of target regions with ONE HSP: $blastingStats{'targets_with_one_hsp'} ("
        . ((($blastingStats{'targets_with_one_hsp'})/($blastIO->result_count))*100) . "% of total)\n";
    print $STATSFH "Number of targets with at least 98% of their length within a contig: $blastingStats{'targets_nested_within_contig'} ("
        . ((($blastingStats{'targets_nested_within_contig'})/($blastIO->result_count))*100) . "% of total)\n";
    
    close($STATSFH) or die "Can't close blast stats summary file handle.\n";
    
    unlink "$velvet_dir_name/Graph2";
    unlink "$velvet_dir_name/LastGraph";
    unlink "$velvet_dir_name/PreGraph";
    unlink "$velvet_dir_name/Sequences";
    unlink "$velvet_dir_name/Roadmaps";
    
    # The maximum number of targets (ymax in the graph) should be $blastIO->result_count();
    # Once done processing, there should be ($toK-$fromK)/2 values along the x-axis,
    # which are elements of the arrays below
    push (@targets_with_hits_for_R, $blastingStats{'targets_with_hits'});
    push (@targets_with_one_hit_for_R, $blastingStats{'targets_with_one_hit'});
    push (@targets_with_one_HSP_for_R, $blastingStats{'targets_with_one_hsp'});
    push (@targets_nested_in_contig_for_R, $blastingStats{'targets_nested_within_contig'});
}
unlink $singlesFile;
unlink $R1_postjoin;
unlink $R2_postjoin;


# Create plots using R
my $plotsFileName = $name . "_plots.jpg";

my $R = Statistics::R->new();
$R->start();
$R->run(q`library("Cairo")`);

$R->set('x',\@kmer_choices);
$R->set('yHits', \@targets_with_hits_for_R);
$R->set('yOneHit', \@targets_with_one_hit_for_R);
$R->set('yOneHSP', \@targets_with_one_HSP_for_R);
$R->set('yNested', \@targets_nested_in_contig_for_R);

$R->run(qq`CairoJPEG("$plotsFileName")`);
$R->run(q`par(mfrow=c(2,2))`);
$R->run(q`plot(x,yHits,type="l",lwd=5,col="blue",xlab="kmer value",ylab="Number of targets",main="Targets with matching contigs")`);
$R->run(q`plot(x,yOneHit,type="l",lwd=5,col="blue",xlab="kmer value",ylab="Number of targets",main="Targets with ONE matching contig")`);
$R->run(q`plot(x,yOneHSP,type="l",lwd=5,col="blue",xlab="kmer value",ylab="Number of targets",main="Targets with ONE matching HSP")`);
$R->run(q`plot(x,yNested,type="l",lwd=5,col="blue",xlab="kmer value",ylab="Number of targets",main="Targets WITHIN a contig")`);
$R->run(q`dev.off()`);

$R->stop();







#Documentation
__END__

=head1 NAME

velvetTarget.pl

=head1 SYNOPSIS 

perl velvetTarget.pl --R1 <file> --R2 <file> --name <string> --from <int> \
--to <int> --adapters <file> --out <file> --probes <file>

 Options:
   -name=s          Name to be used as prefix for directories, etc...
   -R1=s            R1 reads file
   -R2=s            R2 reads file
   -adapters=s      Fasta file containing the adapters to trim
   -from=i          Starting kmer value (default 19)
   -to=i            Ending kmer value (default 201)
   -out=s           Log file name
   -probes=s        A fasta file containing the probes or target regions
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
