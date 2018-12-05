### Functions for Simulating RILS

library('data.table')
library('dplyr')
library('tibble')
setwd('~/Documents/PBGG/MAGIC/')
ogutmap=fread('data_files/ogutmap_v4_ordered.csv',data.table=F)

set.seed(11)
options(scipen=99)

chrom_init <- function(c=10,h1_breakpoints=NULL,h1_donors=c('A'),h2_breakpoints=NULL,h2_donors=c('A')){
  chrom_end=max(ogutmap[ogutmap$chr==c,]$pos)
  if(is.null(h1_breakpoints) & is.null(h2_breakpoints)){
    h1_breakpoints=c(chrom_end)
    h2_breakpoints=c(chrom_end)
  }
  ind = list(chr=c,h1=list(breakpoints=h1_breakpoints,donors=h1_donors),h2=list(breakpoints=h2_breakpoints,donors=h2_donors))
  return(ind)
}

make_f1 <- function(p1,p2,c=10){
  if(runif(1)>=0.5){
    h1_donors=p1$h1$donors
    h2_donors=p2$h1$donors
    h1_breakpoints=p1$h1$breakpoints
    h2_breakpoints=p2$h1$breakpoints
  }
  else{
    h1_donors=p2$h1$donors
    h2_donors=p1$h1$donors
    h1_breakpoints=p2$h1$breakpoints
    h2_breakpoints=p1$h1$breakpoints
  }
  f1=chrom_init(c=c,h1_breakpoints = h1_breakpoints,h1_donors=h1_donors,h2_breakpoints=h2_breakpoints,h2_donors=h2_donors)
  return(f1)
}

no_breaks <- function(donor){
    return(length(donor$breakpoints)==1)
}


interference <- function(xo_pos){
  if(xo_pos[2]-xo_pos[1] < 10e6){
    choice=sample(c(1,2))
    return(xo_pos[choice])
  }
  else{
    return(xo_pos)
  }
}

crossover<-function(parent,recomb,c){
  ### parents: Dataframe with two columns, A pair of parents to use to generate
  ### the F1
  ### c: The chromosome being simulated
  
  ### Returns: a dataframe of n lines that are F1s of the 2 parents
  interp=approxfun(recomb$cumprob,recomb$pos,yleft=min(recomb$pos),yright=max(recomb$pos))
  xo = sample(c(1,2),1,p=c(0.75,0.25))
  draw=runif(xo)
  xo_pos=round(interp(draw),0)
  xo_pos=sort(xo_pos)
  if(length(xo_pos)==2){
    xo_pos=interference(xo_pos)
  }
  if(runif(1)>=0.5){
    donor1=parent$h1
    donor2=parent$h2
    current='h1'
    }
  else{
    donor1=parent$h2
    donor2=parent$h1
    current='h2'
    }
  start=1
  for(p in xo_pos){
    if(no_breaks(donor1) & no_breaks(donor2)){
      recomb_pos=c(p,donor1$breakpoints)
      recomb_donor=c(donor1$donors,donor2$donors)
    }
    else if(p < donor1$breakpoints[1] & p < donor2$breakpoints[1]){
      recomb_pos=c(p,donor2$breakpoints)
      recomb_donor=c(donor1$donors[1],donor2$donors)
    }
    else if(p < donor1$breakpoints[1]){
      d2_upper=donor2$breakpoints[donor2$breakpoints-p>0]
      d2_index=which(min(d2_upper)==d2_upper)
      
      recomb_pos=c(p,donor2$breakpoints[d2_index:length(donor2$breakpoints)])
      recomb_donor=c(donor1$donors[1],donor2$donors[d2_index:length(donor2$donors)])
    }
    else if(no_breaks(donor1)==T & no_breaks(donor2)==F){
      d2_upper=donor2$breakpoints[donor2$breakpoints-p>0]
      d2_index=which(min(d2_upper)==d2_upper)
      recomb_pos=c(p,donor2$breakpoints[d2_index:length(donor2$breakpoints)])
      recomb_donor=c(donor1$donors,donor2$donors[d2_index:length(donor2$donors)])
    }
    else if(no_breaks(donor1)==F & no_breaks(donor2)==T){
      d1_lower=donor1$breakpoints[p-donor1$breakpoints>0]
      if (length(d1_lower)==0){
        d1_index=1
        recomb_pos=c(donor1$breakpoints[1:d1_index],p,donor2$breakpoints)
        recomb_donor=c(donor1$donors[1:d1_index],donor2$donors)
      }
      else{
        d1_index=which(max(d1_lower)==d1_lower)
        recomb_pos=c(donor1$breakpoints[1:d1_index],p,donor2$breakpoints)
        d1_index=d1_index+1
        recomb_donor=c(donor1$donors[1:d1_index],donor2$donors)
      }
    }
    else{
      d1_lower=donor1$breakpoints[p-donor1$breakpoints>0]
      d2_upper=donor2$breakpoints[donor2$breakpoints-p>0]
      d2_index=which(min(d2_upper)==d2_upper)
      if (length(d1_lower)==0){
        d1_index=1
        recomb_pos=c(donor1$breakpoints[1:d1_index],p,donor2$breakpoints[d2_index:length(donor2$breakpoints)])
        recomb_donor=c(donor1$donors[1:d1_index],donor2$donors[d2_index:length(donor2$donors)])
      }
      else{
        d1_index=which(max(d1_lower)==d1_lower)
        recomb_pos=c(donor1$breakpoints[1:d1_index],p,donor2$breakpoints[d2_index:length(donor2$breakpoints)])
        d1_index=d1_index+1
        recomb_donor=c(donor1$donors[1:d1_index],donor2$donors[d2_index:length(donor2$donors)])
      }
    }
    if(current=='h1'){
      donor1=list(breakpoints=recomb_pos,donors=recomb_donor)
      donor2=parent$h1
      current='h2'
    }
    else{
      donor1=list(breakpoints=recomb_pos,donors=recomb_donor)
      donor2=parent$h2
      current='h2'
    }
  }
  return(list(breakpoints=recomb_pos,donors=recomb_donor))
}


get_gametes<-function(parent,recomb,c){
  event=sample(c(1,2,3),1)
  if(event==1){
    result=parent$h1
  }
  else if(event==2){
    result=parent$h2
  }
  else{
    result=crossover(parent,recomb,c)
  }
  return(result)
}

offspring<-function(p1,p2,chroms){
  offspring=list()
  for(c in chroms){
    recomb=ogutmap[ogutmap$chr==c,]
    h1=get_gametes(p1[[c]],recomb,c)
    h2=get_gametes(p2[[c]],recomb,c)
    offspring[[c]]=list(chr=c,h1=h1,h2=h2)
  }
  return(offspring)
}


#Generate two test diploid parents A and B
#Produce an F1 from Aand B

P1=list()
for(x in seq(1,10)){P1[[x]]=chrom_init(x,h1_donor='A',h2_donor='A')}
P2=list()
for(x in seq(1,10)){P2[[x]]=chrom_init(x,h1_donor='B',h2_donor='B')}

f1=list()
for(x in seq(1,10)){f1[[x]]=make_f1(P1[[x]],P2[[x]],x)}

f2=offspring(f1,f1,chroms=seq(1,10))

make_ril <- function(f2,ngen,c){
  n=ngen
  if(n==0){
    return(f2)
  }
  else{
    ril=offspring(f2,f2,c)
    make_ril(ril,ngen=n-1,c)
  }
}

ril=make_ril(f2,6,c=seq(1,10))

nam_founders=c('B97','CML228','CML277','CML333','CML69','IL14H','KI3','M162W','MO18W','NC350','OH43','P39','TZI8',
               'CML103','CML247','CML322','CML52','HP301','KI11','KY21','M37W','MS71','NC358','OH7B','TX303')

make_nam <- function(parents,ngen,c){
  #initiate founders
  #make B73
  nam=list()
  B73=list()
  for(x in seq(1,10)){B73[[x]]=chrom_init(x,h1_donor='B73',h2_donor='B73')}
  founders=list()
  for(f in parents){
    parent=list()
    for(x in seq(1,10)){parent[[x]]=chrom_init(x,h1_donor=f,h2_donor=f)}
    
    f1=list()
    for(x in seq(1,10)){f1[[x]]=make_f1(parent[[x]],B73[[x]],x)}
    f6_pop=list(parent=f)
    for(i in seq(1,200)){
      for( x in seq(1,10)){
        f6_pop[[i]]=list()
        f6_pop[[i]][[x]]=make_ril(f1,ngen,c)
      }
    }
    nam[[f]]=f6_pop
  }
  return(nam)
}

nam = make_nam(nam_founders,6,seq(1,10))
