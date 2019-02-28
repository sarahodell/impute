#!/env/usr/bin Rscript

library('data.table')
library('dplyr')
library('tibble')

args=commandArgs(trailingOnly=TRUE)
c=as.character(args[1])


founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra","FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")

binfinal=fread(sprintf('wgs_founders/Biogemma_WGS_all_alleles_final_chr%s.txt',c),data.table=F)

#Turn all NA values to zero for calculating allele probs.
freq=rowMeans(binfinal[,5:20],na.rm=T)
for(i in 5:20){
    ind=which(is.na(binfinal[,i]))
    binfinal[ind,i]=freq[ind]
}

options(scipen=999)

#for each chromosome
#read in files
pr=readRDS(sprintf('data/bg%s_genoprobs.rds',c))
#Read in physical map
pmap=fread(sprintf('qtl2_startfiles/Biogemma_pmap_c%s.csv',c),data.table=F)

wgslen=dim(binfinal)[1]
sub=as.matrix(binfinal[,c(founders)])
pos=binfinal$pos
markers=binfinal$marker
samples=names(pr[[1]][,1,1])
alt1=binfinal$alt1
ref=binfinal$ref
remove(binfinal)
all_probs=c()
for( i in seq(1,344)){
    print(i)
    ind=pr[[1]][i,,]
    ind=t(ind)
    ind=as.data.frame(ind)
    ind=rownames_to_column(ind,"marker")
    names(ind)=c("marker",founders)
    ind=merge(ind,pmap,by.x="marker",by.y="marker")
    ind=ind[order(ind$pos),]
    rownames(ind)=seq(1,nrow(ind))
    ind=ind[,c("marker","chr","pos",founders)]
    prlen=dim(ind)[1]
    f_probs=c()
    for(f in seq(1,16)){
        founder=founders[f]
        f_ind=ind[,c('pos',founder)]
        f_interp=approxfun(f_ind$pos,f_ind[,c(founder)],method="linear",yleft=unlist(unname(f_ind[1,c(founder)])),yright=unlist(unname(f_ind[prlen,c(founder)])))
        f_probs=rbind(f_probs,f_interp(pos))
    }
    f_probs=t(as.matrix(f_probs))
    allele_probs=rowSums(sub*f_probs)
    all_probs=rbind(all_probs,allele_probs)
}
all_probs_t=t(all_probs)
remove(all_probs)
all_probs_t=as.data.frame(all_probs_t)
rownames(all_probs_t)=markers
names(all_probs_t)=samples
all_probs_t=rownames_to_column(all_probs_t,"marker")
all_probs_t$alt1=alt1
all_probs_t$ref=ref
all_probs_t=all_probs_t[,c('marker','alt1','ref',samples)]
print("Writing allele probs to file")
fwrite(all_probs_t,sprintf('data/bg%s_wgs_alleleprobs.txt',c),row.names=F,sep='\t',quote=F)