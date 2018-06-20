#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys


def arg_parse():
    parser=argparse.ArgumentParser(description="""full_sim.py takes simulates a MAGIC or RIL population""")
    parser.add_argument("-f",type=str,help="""file containing names of parents to build lines from (one per line)""")
    parser.add_argument("-n",type=int,help="""Number of lines to procuce (int)""")
    parser.add_argument("-t",type=str,help="""Type of lines to simulate. Options are 'magic' or 'ril'""")
    args = parser.parse_args()
    return args


def chrom_info(chro_num):
    """chro_num: chromosome number (int)
    Outputs:
    end: (int) the approximate end of the chromosome, in Mb
    """
    cent = pd.read_table('B73v4centromeres.txt',sep='\t')
    cdf = cent[cent['chr']==chro_num]
    end = int(cdf['v3chr.end'].values[0])
    centstart = int(cdf['v3start'].values[0])
    centend= int(cdf['v3end'].values[0])
    centlen = int(round(cdf['v3size'].values[0]))
    if end-centend > centstart:
        longarm = end-centend
        shortarm = centstart
    else:
        longarm = centstart
        shortarm = end-centend
    diff = longarm-shortarm
    a = [i**2 for i in list(reversed(range(0,centstart+1)))]
    c = [0.0 for i in range(centlen)]
    b = [i**2 for i in range(longarm)]
    xo_prob = a+c+b
    if len(xo_prob)>= end:
        xo_prob=xo_prob[:end]
    xo_prob = [float(j)/sum(xo_prob) for j in xo_prob]
    return xo_prob,end


def chrom_sim(founders,pnum,c=10):
    """
    pnum: (int) Number of parents to start with. Parents are randomly selected
    from the 26 NAM founders
    c: chromosome number (1..10 for maize)
    Outputs a 2-D list of length pnum. Internal lists contain strings 
    assigning parental donor for 1Mb blocks (i.e. list[0]= ['B73','B73',...])"""
    pop=[]
    if len(founders)==pnum:
        parents=founders
    else:
        parents = np.random.choice(founders,pnum,replace=False)
    for i in parents:
        xo_prob,end=chrom_info(c)
        p = []
        for j in range(end+1):
            p.append(i)
        pop.append(p)
    return pop



def crossover(n,parents,c=10):
    """
    Simulates a crossover event of a chromosome, returns the 'f1'
    Input:
    n: (int) number of f1 samples to produce
    parents: 2-D list of length 2, output of chrom_sim() (i.e. chromsim(16)[:2])
    c: chromosome number 
    
    Output:
    If n==1, 1-D list with f1
    If n>1, 2-D list of f1s
    """
    rils = []
    for i in range(n):
        prior,end=chrom_info(c)
        site = [s for s in range(end)]
        f1=[]
        #randomly choose a parent to start with
        if np.random.random_sample() >= 0.5:
            donor = parents[0]
        else:
            donor = parents[1]
        xo = np.random.choice([1,2],p=[0.6,0.4])
        draw = np.random.choice(site,size=xo,p=prior)
        if xo == 2:
            while abs(draw[1] - draw[0]) < 40:
                draw = [draw[0],np.random.choice(site,p=prior)]
        draw = sorted(draw)
        #iterate through and take sections from each parent
        start = 0
        for d in draw:
            f1+=donor[start:d]
            start=d
            if donor==parents[0]:
                donor=parents[1]
            elif donor==parents[1]:
                donor=parents[0]
        f1+=donor[start:]
        rils.append(f1)
    if n==1:
        return f1
    else:
        return rils



def make_magic(parents,c=10,n=1):
    """ Simulates MAGIC lines
    Input:
    parents: (list)2-D list of parent chromosomes (length must be even
    output of chrom_sim
    c: chromosome number (Default: 10)
    n: number of f1 samples to make per parent cross (Default: 1)
    
    Output:
    2-D list of length(n), simulated chromosomes
    """
    if len(parents)==2:
        return crossover(n,parents)
    rounds = []
    for i in range(0,len(parents),2):
        rounds.append(crossover(n,parents[i:i+2]))
    return make_magic(rounds)


def locations(r,c=10):
    """Identifies chromosome breakpoints and returns a list of format:
    [[chr,start,end,donor],...]
    Input: (list) simulated chromosome (i.e. f[0])
    """
    locs=[]
    last = r[0]
    counter = 0
    for i in range(len(r)):
        if r[i] != last:
            start = i
            locs.append([c,counter*1e6,(start*1e6)-1,last])
            counter=start
            last = r[i]
    locs.append([c,start*1e6,len(r)*1e6,last])
    return locs


def make_outfile(ril,out):
    """ Writes out file in format: 
    out: name of outfile (str)
    """
    txt = 'sample\tchr\tstart\tend\tdonor1\tdonor2\n'
    count=1
    for j in ril:
        locs = locations(j)
        for l in locs:
            txt+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('M'+str(count),l[0],int(l[1]),int(l[2]),l[3],l[3])
        count+=1
    with open(out,'w') as outfile:
        outfile.write(txt)



        
if __name__ == "__main__":
    pass
