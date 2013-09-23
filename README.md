velvetTarget
============

Choosing the right kmer value is critical for de novo assembly of contigs from next generation sequencing, but this is difficult. This task is easier for genome reconstruction, as the goal is usually to assemble the longest contigs possible (VelvetOptimiser.pl is great for this). For target enrichment experiments, however, the goal is rather to recover as many of the targets as completely as possible.

This script attempts to find the k value that achieves this. Instead of basing the choice of k off of summary statistics of the assembled contigs themselves, it rather blasts the baits or target regions against the assembly, and returns the k value and assembly that recovers the most targets the most completely.

Assumptions:
    *You have paired end HiSeq or MiSeq data
    *You have three data files--R1, R2, and singletons/joined reads
    *You have compiled velvet with a max kmer value at least as high as you want to explore (program default goes up to 201)