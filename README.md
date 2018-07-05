# Imputation Validation Using Simulated Crosses

Generate simulated recombinant lines based off of the NAM maize population 
(5 generations of selfing) and MAGIC data from Biogemma. Writes out into files for other building simulated vcf files or for visualization in the plots below.

Code to visualize parent donors of RIL chromosomes in an ideogram. This is modified from
 code by Ryan Dale (https://gist.github.com/daler/c98fc410282d7570efc3)

SimulateRILs: Code for generating and visualizing simulated RILs

FILLINStart: Code for setting up TASSEL, using FILLIN to build haplotypes and impute GBS v2.7 data for the NAM RILs using the NAM founders in HapMap v3.21

FILLINImputeSims: Code for generating VCF files for simulated RILs, imputing these RILs with some proportion of data missing using FILLIN and the RIL parent HapMap v3.21 data

ValidateFILLINPipeline: Full pipeline for creating, imputing, and assessing accuracy of FILLIN for MAGIC lines 

data_files: directory contains files that are

sim_files: directory contains files that have been simulated, ordered by date generated

pscripts: directory contains Python scripts that are called in the Jupyter notebooks above


qtl2: directory contains scripts and data files from attempts to impute with the software R/qtl2 and for formatting the input files for this software.



