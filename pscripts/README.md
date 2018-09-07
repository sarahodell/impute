### pscripts directory

- build_simvcf.py : script takes in a generated files of crossover locations and parental donors and constructs simulated vcf files from donor files

- full_sim.py : module or script for full pipeline of magic line simulation

- intersect.py : Takes in the actual simulated file and the predicted file from FILLIN and calculates proportion of chromosome correctly assigned to the right parent

- marker_generator.py : Randomly selects a set number of markers from a vcf file to use for constructing vcf files with build_simvcf.py

- PCA_600K.Rmd : Markdown of Principle Component Analysis of the 600K Genotype data for the founder and DH lines using SNPRelate in R