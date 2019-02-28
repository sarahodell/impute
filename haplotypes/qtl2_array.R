#!/usr/bin/env Rscript
###Run R/qtl2 on 344 16-way MAGIC DH lines 

library('qtl2')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])

bg<-read_cross2(sprintf('Biogemma_c%s.json',c))

pr <- calc_genoprob(bg,error_prob=0.002,cores=4)
print(dim(pr[[1]]))
saveRDS(pr,sprintf("bg%s_genoprobs_010319.rds",c))