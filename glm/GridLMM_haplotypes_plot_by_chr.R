#!/usr/bin/env Rscript

#### Run GridLMM on a phenotype for and environment across all h haplotype groups for one chromosome

args=commandArgs(trailingOnly=T)
pheno=as.character(args[[1]])
env=as.character(args[[2]])
#chr=as.character(args[[3]])
#cores=as.numeric(args[[4]])

#date=format(Sys.time(),'%m%d%y')

library('ggplot2')
library('data.table')
library('dplyr')

all_chroms=c()

for(i in 1:10){
   pmap=fread(sprintf('../qtl2_startfiles/Biogemma_pmap_c%.0f.csv',i),data.table=F)
   ID=c()
   pos=c()
   bin=c()
   pvalues=c()
   hapgrp=c()
   for(h in 2:16){
      mod=readRDS(sprintf('models/Biogemma_chr%.0f_haplogrp%.0f_%s_x_%s.rds',i,h,pheno,env))
      p=match(mod$X_ID,pmap$marker)
      phy_pos=pmap[p,]$pos
      ID=c(ID,mod$X_ID)
      pos=c(pos,phy_pos)
      pvalues=c(pvalues,mod$p_value_ML)
      hapgrp=c(hapgrp,rep(h,length(mod$p_value_ML)))
      bin=c(bin,seq(1,length(mod$p_value_ML)))
   }
   gwas=data.frame(chr=i,hapgrp=hapgrp,bin=bin,ID=ID,pos=pos,pvalues=pvalues,stringsAsFactors=F)
   all_chroms=rbind(all_chroms,gwas)
}
all_chroms$chr=as.factor(all_chroms$chr)
all_chroms$bin=as.factor(all_chroms$bin)
all_chroms$hapgrp=as.factor(all_chroms$hapgrp)
all_chroms$log10p = -log10(all_chroms$pvalues)

fwrite(all_chroms,sprintf('sig_tables/%s_x_%s.txt',pheno,env),quote=F,row.names=F,sep='\t')
size=dim(all_chroms)[1]
#cutoff=3.521018
cutoff=6.037128

all_chroms$sig = all_chroms$log10p>= cutoff

png(sprintf('images/%s_x_%s_manhattan_sig.png',pheno,env),width=960,height=680)
theme_set(theme_classic())
theme_update(text=element_text(family="Helvetica"))
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(size=26),axis.title=element_text(size=14,face="bold"))
theme_update(panel.background=element_blank())

print(ggplot(all_chroms,aes(x=pos/1e6,y=log10p)) + geom_point(aes(color=sig)) + scale_color_manual(breaks=all_chroms$sig,values=c("FALSE"="black","TRUE"="blue")) + facet_grid(.~chr,scales="free") + ggtitle(sprintf("%s in %s Using Haplotype Probabilities",pheno,env)) + xlab("Chromosome (Mb)") + ylab("-log10(P-Value)") + guides(color=F))
dev.off()



