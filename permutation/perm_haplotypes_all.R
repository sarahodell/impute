#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
chr=as.character(args[[2]])
h=as.numeric(args[[3]])
cores=as.numeric(args[[4]])
reps=as.numeric(args[[5]])

library('GridLMM')
library('data.table')
library('dplyr')
library('lme4')
library('parallel')
library('MASS')

all_reps=c()

# Read in Kinship Matrix
K=fread(sprintf('../K_matrices/K_matrix_chr%s.txt',chr),data.table=F)
rownames(K)=K[,1]
rownames(K)=gsub("-",".",rownames(K))
K=as.matrix(K[,-1])
colnames(K)=rownames(K)

# Read in phenotypes
# Grab the phenotype of interest and drop the genotypes not in the K matrix
phenotypes=fread('../phenotypes.csv',data.table=F)
phenotypes=phenotypes[,c('Genotype_code','Loc.Year.Treat',pheno)]
phenotypes$Genotype_code=gsub('-','.',phenotypes$Genotype_code)
phenotypes=phenotypes[phenotypes$Genotype_code %in% rownames(K),]

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Loc.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)

m1=lmer(y~Loc.Year.Treat + (1|ID2),data)
data_blup = as.data.frame(ranef(m1)$ID2)
data_blup$ID = rownames(data_blup)
data_blup$y=data_blup$`(Intercept)`
data_blup=data_blup[,c('ID','y')]

null_model = GridLMM_ML(y~1 + (1|ID),data_blup,relmat=list(ID=K),ML=T,REML=F,verbose=F)

# Read in the haplotype group probabilities
# Filter genotypes that are not in the K matrix
X_list=readRDS(sprintf('../../haplotype_probs/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))
X_list_full=lapply(X_list,function(x) x[data_blup$ID,])

n_reps=seq(1,reps)

randomized_gwas <- function(rep) {

   ### Randomly reassign genotypes to individuals
   len=dim(X_list_full[[1]])[1]
   draw=sample(len,len,replace=F)
   
   X_list_reordered=lapply(X_list_full,function(x) x[draw,])
   for(x in seq(1,h)){
       dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_full[[1]])[[1]]
   }

   # Run GridLMM

   h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
   names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
   h2_start
   V_setup=null_model$setup

   Y=as.matrix(data_blup$y)
   X_cov=null_model$lmod$X
   X_list_null=NULL

   gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
   gwas=gwas[!is.na(gwas$p_value_ML),]
   tmp=data.frame(chr=chr,hapgrp=h,replicate=rep,pval=min(gwas$p_value_ML))
}

print(system.time({
results=mclapply(n_reps,randomized_gwas,mc.cores=cores)
}))

saveRDS(results,sprintf('test_models/chr%s_haplogroup%.0f_%s_x_ALL_%.0frep_max_pvalues.rds',chr,h,pheno,reps))
