Folder contains files for contructing input files for R/qtl2

foundergeno.py - Creates a founder genotypes csv file from an inputed csv file
	       - input csv file can be created via


``` bash

module load bcftools
header='id,chr,pos,[...samples]

echo $header > tmp.csv
bcftools query -f '%ID,%CHROM,%POS[,%GT]\n' founder_genotypes.vcf.gz >> tmp.csv

python foundergeno.py tmp.csv outfile.csv
```


genofile2.py - Creates a genotype csv file from input vcf files for samples of interest

pmap.py - Creates a physical map csv file from input vcf file with markers and and physical positions being used for qtl2
