---
title: "FILLINImputeSims"
author: "Sarah Odell"
date: "5/14/2018"
output: html_document
---

#Impute Simulated RILs Using FILLIN

###Getting Started
Software Requirements:
Tassel 5.0
bcftools/1.2
tabix/0.2.6

Summary :
Create vcf files for simulated RILs (B73 x Oh43 in this example) and drop data. 
Impute the RILs with some proportion of data missing and validate how reliably the proper
 parent is imputed onto the RIL.

For information on how RILs were generated, see SimulateRILs. For information on how to
 set up TASSEL see FILLINStart.

The file 'RILS_10x10.txt' was generated using the code in SimulateRILs. It contains 
breakpoint information for 10 simulated RILs for chromosome 10. Using this file, we will
 generate a vcf file for the simulated RILs.

###Preparing VCF Files
The file c10_hmp31_edit_founders.vcf.gz should be indexed already. In this example, the RIL
parents were B73 and Oh43, so we will make individual vcf files for both of them
to make extracting variant information from both of them easier

```{bash extractparents}
for i in B73 Oh43; do
    bcftools view -Oz -s $i  c10_hmp31_edit_founders.vcf.gz > ${i}_c10_hmp321.vcf.gz
    echo $i>>parents.txt 
done

### Then combine into one file (will be used later for building haplotypes)
bcftools view -Oz -S parents.txt c10_hmp31_edit_founders.vcf.gz > RILparents_c10_hmp321.vcf.gz 
```

###Building RIL VCF Files
The script build_ril.py takes the following arguments:

```{bash buildril}
python build_ril.py --help
```

To generate a RIL file without dropping data, just run:

```{bash fullril}
python build_ril.py RILS_10x10.txt B73xOh43_10samples.vcf
```

This will create a vcf file for each of the 10 RIL samples. The following code reheaders
 the files with the proper sample names and merges them all into one file. I am working on
  a less sloppy way to do this.
  
```{bash rilmerge}
for i in 1 2 3 4 5 6 7 8 9 10; do
    echo R${i} > R${i}.txt
    bcftools reheader -s R${i}.txt R${i}_B73xOh43_10samples.vcf > R${i}_B73xOh43_10samples_edit.vcf
    bgzip R${i}_B73xOh43_10samples_edit.vcf
    tabix -p vcf R${i}_B73xOh43_10samples_edit.vcf.gz
    echo R${i}_B73xOh43_10samples_edit.vcf.gz >> ril_names.txt
done

bcftools merge -l ril_names.txt -m all -Oz > B73xOh43_RILSimsAll_chr10.vcf.gz
```

###Building RIL VCF Files with Missing Data

To generate a RIL file with 40% of data missing, use the options --drop True --droprate 0.4:

```{bash dropril}
python build_ril.py RILS_10x10.txt B73xOh43_10samples_drop40.vcf --drop True --droprate 0.4
```

This will create a vcf file for each of the 10 RIL samples. The following code reheaders
 the files with the proper sample names and merges them all into one file.
 
```{bash dropmerge}
for i in 1 2 3 4 5 6 7 8 9 10; do
    echo R${i} > R${i}.txt
    awk -v OFS="\t" '$1=$1' R${i}_B73xOh43_10samples_drop40_reheader.vcf | bcftools reheader -s R${i}.txt | bgzip -c > R${i}_B73xOh43_10samples_drop40.vcf.gz
    tabix -p vcf R${i}_B73xOh43_10samples_drop40.vcf.gz
    echo R${i}_B73xOh43_10samples_drop40.vcf.gz >> ril_names.txt
done

bcftools merge -l ril_names.txt -m all -Oz > B73xOh43_RILSims_chr10_drop40.vcf.gz
```

### Imputing onto the dropped files
Use the TASSEL FILLINFindHaplotypesPlugin to build haplotypes from the RIL parents, Oh43 and B73

```{bash findhaps}
run_pipeline.pl -Xmx32G -FILLINFindHaplotypesPlugin -hmp RILparents_c10_hmp321.vcf.gz -o B73xOh43Haps/ -mxErr 0 -minTaxa 1 -extOut true -endplugin
```

Use the TASSEL FILLINImputationPlugin to impute onto the simulated RIL files with missing data,
and get a projection alignment file

```{bash impute}
run_pipeline.pl -Xmx32G -FILLINImputationPlugin -hmp B73xOh43_RILSims_chr10_drop40.vcf.gz -d RILparents_c10_hmp321.vcf.gz -o B73xOh43_RILSims_chr10_drop40_FILLINImputed.vcf.gz -ProjA true -endPlugin

```

This will impute the parent variants onto the simulated RILS.
Now, we must validate that the proper parent imputed in the proper chromosomal locations








