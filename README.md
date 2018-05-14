# Simulate RILs

Generate recombinant inbred lines (RILs) based off of the NAM maize population 
(5 generations of selfing). Writes out into files for other building RIL vcf files or for 
visualization in the plots below.

Code to visualize parent donors of RIL chromosomes in an ideogram. This is modified from
 code by Ryan Dale (https://gist.github.com/daler/c98fc410282d7570efc3)

SimulateRILs: Code for generating and visualizing simulated RILs
B73v3centromres.txt: File used for generating RILs. Chromosome length and centromere position information based off of the B73v4 maize reference genome (https://www.nature.com/articles/nature22971).
R1.txt: Ideogram file used for plotting the RIL ideogram (randomly generated, so yours will be different)