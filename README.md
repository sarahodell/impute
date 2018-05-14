# Simulate RILs

Generate recombinant inbred lines (RILs) based off of the NAM maize population 
(5 generations of selfing). Writes out into files for other building RIL vcf files or for 
visualization in the plots below.

Code to visualize parent donors of RIL chromosomes in an ideogram. This is modified from
 code by Ryan Dale (https://gist.github.com/daler/c98fc410282d7570efc3)

SimulateRILs: Code for generating and visualizing simulated RILs

FILLINStart: Code for setting up TASSEL, using FILLIN to build haplotypes and impute GBS v2.7 data for the NAM RILs using the NAM founders in HapMap v3.21

FILLINImputeSims: Code for generating VCF files for simulated RILs, imputing these RILs with some proportion of data missing using FILLIN and the RIL parent HapMap v3.21 data

B73v4centromres.txt: File used for generating RILs. Chromosome length and centromere position information based off of the B73v4 maize reference genome (https://www.nature.com/articles/nature22971).

R1.txt: Ideogram file used for plotting the RIL ideogram (randomly generated, so yours will be different)

RILS_10x10.txt File containing breakpoint and parent donor information for chromosome 10 for 10 simulated RILs.
