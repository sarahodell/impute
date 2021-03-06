#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
chr=as.character(args[[3]])
h=as.numeric(args[[4]])
cores=as.numeric(args[[5]])
reps=as.numeric(args[[6]])

library('GridLMM')
library('data.table')
library('dplyr')
library('parallel')
library('MASS')

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

X_list=readRDS(sprintf('../../haplotype_probs/bg%s_filtered_haplogroup%.0f_probs.rds',chr,h))

data=data.frame(ID=phenotypes$Genotype_code,ID2=phenotypes$Genotype_code,Loc.Year.Treat=phenotypes$Lo\
c.Year.Treat,y=phenotypes[,c(pheno)],stringsAsFactors=F)
data=data[data$Loc.Year.Treat==env,]
data=data[!is.na(data$y),]
data$y=data$y - mean(data$y)

null_model = GridLMM_ML(y~1+(1|ID),data,relmat=list(ID=K),ML=T,REML=F,verbose=F)

# Grab only genotypes in data
X_list_full=lapply(X_list,function(x) x[data$ID,])

n_reps=seq(1,reps)

randomized_gwas<-function(rep){
   len=dim(X_list_full[[1]])[1]

   # Run GridLMM

   # randomize the order of the genotypes
   draw=sample(len,len,replace=F)
   X_list_reordered=lapply(X_list_full,function(x) x[draw,])
   for(x in seq(1,h)){
       dimnames(X_list_reordered[[x]])[[1]]=dimnames(X_list_full[[1]])[[1]]
   }

   h2_start=null_model$results[,grepl('.ML',colnames(null_model$results),fixed=T),drop=FALSE]
   names(h2_start) = sapply(names(h2_start),function(x) strsplit(x,'.',fixed=T)[[1]][1])
   h2_start
   V_setup=null_model$setup

   Y=as.matrix(data$y)
   X_cov=null_model$lmod$X
   X_list_null=NULL

   gwas=run_GridLMM_GWAS(Y,X_cov,X_list_reordered[-1],X_list_null,V_setup=V_setup,h2_start=h2_start,method='ML',mc.cores=cores,verbose=F)
   gwas=gwas[!is.na(gwas$p_value_ML),]
   tmp=data.frame(chr=chr,hapgrp=h,replicate=b,pval=min(gwas$p_value_ML))
}

print(system.time({
results=mclapply(n_reps,randomized_gwas,mc.cores=cores)
}))

saveRDS(results,sprintf('test_models/chr%s_haplogroup%.0f_%s_x_%s_%.0frep_max_pvalues.rds',chr,h,pheno,env,reps))


