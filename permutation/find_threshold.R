#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])

library('data.table')
library('dplyr')
library('ggplot2')

options(scipen=20)

all_reps=c()
for(c in 1:10){
   haps=seq(2,16)
   for(h in haps){
         r=readRDS(sprintf('test_models/chr%.0f_haplogroup%.0f_%s_x_%s_1000rep_max_pvalues.rds',c,h,pheno,env))
	 df=c()
	 df=sapply(seq(1,1000),function(x) rbind(df,unlist(r[[x]])))
	 df=t(df)
	 df=as.data.frame(df)
	 names(df)=c('chr','hapgrp','replicate','pval')
	 #df$chr=c
	 #df$replicate=as.numeric(df$replicate)
	 df=df[!is.na(df$pval),]
         #tmp=data.frame(chr=c,replicate=df$replicate,hapgrp=h,pval=df$pval,stringsAsFactors=F)
         all_reps=rbind(all_reps,df)
   }
}

fwrite(all_reps,sprintf('max_reps/%s_x_%s_rep1000_max_pvalues.txt',pheno,env),quote=F,row.names=F,sep='\t')

minp = all_reps %>% group_by(replicate) %>% summarize(pval=min(pval))
minp=as.data.frame(minp)

threshold=quantile(minp$pval,0.05,lower.tail=T)
print(threshold)
print(-log10(threshold))

png(sprintf('%s_x_%s_perm_1000_pval_dist.png',pheno,env))
print(ggplot(minp,aes(x=pval)) + geom_histogram() + geom_vline(xintercept=threshold))
dev.off()

png(sprintf('%s_x_%s_perm_1000_log10pval_dist.png',pheno,env))
print(ggplot(minp,aes(x=-log10(pval))) + geom_histogram() + geom_vline(xintercept=-log10(threshold)))
dev.off()





