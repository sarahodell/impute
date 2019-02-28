#!/usr/bin/env Rscript

library("qtl2")
library("data.table")
library("ggplot2")
library("reshape2")
library("tibble")
library('igraph')
library('abind')

args=commandArgs(trailingOnly=T)
jsonfile=as.character(args[1])
outfile=as.character(args[2])

print(c)
bg=read_cross2(jsonfile)
ibd=find_ibd_segments(bg$founder_geno,bg$pmap,min_lod=15,error_prob=0.002,cores=4)
fwrite(ibd,outfile,quote=F,row.names=F,sep='\t')


