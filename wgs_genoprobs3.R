library('data.table')
library('dplyr')
library('tibble')
setwd('~/Documents/PBGG/MAGIC/qtl2/Biogemma_071118/')

founders=c("A632_usa","B73_inra","CO255_inra","FV252_inra","OH43_inra","A654_inra",
           "FV2_inra","C103_inra","EP1_inra","D105_inra","W117_inra","B96","DK63","F492","ND245","VA85")



#array=fread('Biogemma_600K_founder_alleles_chr10.txt',data.table=F,header=F)
#names(array)=c('chr','pos','ref','alt1','alt2',founders)
#array=array %>% mutate(marker=paste0('S',chr,'_',pos))


#Read in genotypes for the 16 Founders WGS data (Chromosome 10)
wgs=fread('Biogemma_WGS_founder_alleles_chr10.txt',data.table=F,header=F)
names(wgs)=c('chr','pos','ref','alt1','alt2',founders)
wgs=wgs %>% mutate(marker=paste0('S',chr,'_',pos))
wgs=wgs[,c('marker','chr','pos','ref','alt1','alt2',founders)]

#Find sites where there is more than one alternate allele
multiallelic=wgs[wgs$alt2!='.',]
distinct = sapply(seq(1,nrow(multiallelic)), function(x) n_distinct(unname(unlist(multiallelic[x,7:22])),na.rm=T)<=2)
sum(distinct)

#Testing dropping of discordant sites
multiallelic[c(1,2,5),'A654_inra']=NA
#It works!
multiallelic[c(3,6,7),c('A654_inra')]=NA
sites=c(1,2,5,NA,12)
sites=subset(sites,is.na(sites)==F)
#None of the sites are just 1&2 (all of them have reference allele)
#There are only ~86K markers that are multiallelic (not just 0&1) -> Just drop them

wgs=wgs[wgs$alt2=='.',]
dim(wgs)


#t=match(array$marker,wgs$marker)
#final=subset(wgs,!(rownames(wgs) %in% t))
#final=rbind(final,array)
#final=final[order(final$pos),]
#rownames(final)=seq(1,nrow(final))
#final=final[,c('marker','chr','pos','ref','alt1','alt2',founders)]

#Convert to binary 0 or 1
binfinal=ifelse(wgs[,7:22]=='0/0',0,ifelse(wgs[,7:22]=='1/1',1,ifelse(wgs[,7:22]=='2/2',2,NA)))
binfinal=as.data.frame(binfinal)
names(binfinal)=founders
binfinal$marker=wgs$marker
#Drop rows with more than 12 founders with missing data (dropped 27K sites)
binfinal = binfinal[rowSums(is.na(binfinal))<=12,]


binfinal = binfinal %>% mutate(pos=strsplit(marker,'_')[[1]][2])
binfinal$pos=as.numeric(binfinal$pos)
binfinal=binfinal[,c('marker','pos',founders)]


fwrite(binfinal,'Biogemma_founder_alleles_chr10.txt',row.names = F,quote=F,sep='\t')

pr=readRDS('bg10_geno_probs.rds')





drop_missing <- function(alleles,probs){
  x=alleles[!is.na(alleles)]
  y=probs[!is.na(alleles)]*(1/sum(probs[!is.na(alleles)]))
  return(sum(x*y))
}

x=unlist(unname(binfinal[1,3:18]))
y=f_probs[1,1:16]
start_time <- Sys.time()
drop_missing(x,y)
end_time <- Sys.time()

(end_time - start_time)*dim(binfinal)[1]
#multiple_alleles <- function(alleles,probs){
#  x=ifelse(alleles==2,1,0)
#  return(sum(x*probs))
#}
pmap10=fread('startfiles/Biogemma_pmap_c10.csv',data.table=F)
pmap10$pos=pmap10$pos*1e6

test=pr[[1]][1,,]
ind=t(test)
ind=as.data.frame(ind)
ind<-rownames_to_column(ind,"marker")
names(ind)=c("marker",founders)
ind <- merge(ind,pmap10,by.x='marker',by.y='marker')
ind=ind[order(ind$pos),]
rownames(ind)=seq(1,nrow(ind))
ind = ind %>% mutate(marker2=paste0('S',chr,'_',pos))
ind=ind[,c('marker','marker2','chr','pos',founders)]

len=dim(test)[2]
f_probs=c()
for(f in seq(1,16)){
  founder=founders[f]
  f_ind=ind[,c('pos',founder)]
  f_interp=approxfun(f_ind$pos,f_ind[,c(founder)],method='linear',yleft=f_ind[,c(founder)][1],yright=f_ind[,c(founder)][len])
  f_probs=rbind(f_probs,f_interp(binfinal$pos))
}
wgslen=dim(binfinal)[1]
#Probability of allele 1
all_probs=c()
allele_probs=sapply(seq(1,wgsdlen),function(x) drop_missing(unname(unlist(binfinal[x,2:17])),f_probs[1:16,x]))
all_probs=rbind(all_probs,allele_probs)

names(all_probs)=binfinal$marker

start_time <- Sys.time()
allele_probs=sapply(seq(1,100),function(x) drop_missing(unname(unlist(binfinal[x,2:17])),f_probs[1:16,x]))
end_time <- Sys.time()

nas=is.na(binfinal)==T
binfinal[nas]=0
all_probs=c()
for( i in seq(1,10)){
  ind=pr[[1]][i,,]
  ind=t(ind)
  ind=as.data.frame(ind)
  ind<-rownames_to_column(ind,"marker")
  names(ind)=c("marker",founders)
  ind <- merge(ind,pmap10,by.x='marker',by.y='marker')
  ind=ind[order(ind$pos),]
  rownames(ind)=seq(1,nrow(ind))
  ind = ind %>% mutate(marker2=paste0('S',chr,'_',pos))
  ind=ind[,c('marker','marker2','chr','pos',founders)]
  
  len=dim(test)[2]
  f_probs=c()
  for(f in seq(1,16)){
    founder=founders[f]
    f_ind=ind[,c('pos',founder)]
    f_interp=approxfun(f_ind$pos,f_ind[,c(founder)],method='linear',yleft=f_ind[,c(founder)][1],yright=f_ind[,c(founder)][len])
    f_probs=rbind(f_probs,f_interp(binfinal$pos))
  }
  
  #Probability of allele 1
  n=dim(f_probs)[2]
  sub=as.matrix(binfinal[,c(founders)])
  f_probs=t(as.matrix(f_probs))
  #allele_probs=sapply(seq(1,n),function(x) sum(unname(unlist(binfinal[x,2:17]))*f_probs[1:16,x]))
  allele_probs=rowSums(sub*f_probs)
  all_probs=rbind(all_probs,allele_probs)
  
}
all_probs=as.data.frame(all_probs)
names(all_probs)=binfinal$marker
rownames(all_probs)=names(pr[[1]][,1,1])

