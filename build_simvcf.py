#!/usr/bin/env python

from subprocess import Popen, PIPE
import pandas as pd
import sys
import argparse
import numpy as np

def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output vcf file""")
    parser.add_argument("markerfile",type=str,help=""""File with list of marker postions""")
    args = parser.parse_args()
    return args


def co_loc(sample,bedfile,c=10):
    """Returns regions of a chromosome received from particular parents in an individual RIL chromosome
    Input:
    breaks: pandas df of format: "start" "donor"
    c: chromosome number (int)

    Output:
    2-D list of form at [[chromsome,start,end,donor],...]
    """
    s = bedfile[bedfile['sample']==sample]
    locs=[]
    parents = s['donor'].unique()
    for index,row in s.iterrows():
        locs.append([row['chr'],int(row['start']),int(row['end']),row['donor']])
    return locs,parents

    
def marker_regions(pbreaks,markerfile,rfile,c=10):
    markers=[]
    with open(markerfile,'r') as infile:
        for line in infile:
            markers.append(line[:-1])
    regions=""
    for m in markers:
        for i in pbreaks:
            chrom=c
            start=i[1]
            end=i[2]
            if int(m) >= start and int(m) <=end:
                regions+='{0}\t{1}\n'.format(chrom,m)
    with open(rfile,'w') as outfile:
    	outfile.write(regions)


def bcftools_view(donorfile,regionsfile=None,header=False):
    if header ==True:
        process = Popen(['bcftools','view','-h',donorfile],stdout=PIPE,stderr=PIPE)
    else:
        process = Popen(['bcftools','view','-H','-R',regionsfile,donorfile],stdout=PIPE,stderr=PIPE)
    stdout,stderr = process.communicate()
    return stdout,stderr


def main():
    args=get_args()
    bedfile = pd.read_table('{0}'.format(args.infile),sep='\t')
    samples = bedfile['sample'].unique()
    header,stderr=bcftools_view(donorfile='../hmp3_founders2/hmp3_founders_final.vcf.gz',header=True)
    print(stderr)
    for sample in samples:
        vcf = header
        breaks,parents = co_loc(sample,bedfile)
        for i in parents:
            pbreaks = [j for j in breaks if j[3]==i]        
            regionsfile='{0}_{1}_regions.txt'.format(i,sample)
            marker_regions(pbreaks=pbreaks,markerfile=args.markerfile,rfile=regionsfile,c=10)
            positions,stderr=bcftools_view(donorfile='../hmp3_founders2/{0}_c10_hmp321_final.vcf.gz'.format(i),regionsfile=regionsfile)
            print(stderr)
	    vcf+=positions
        with open('{0}_{1}'.format(sample,args.outfile),'w') as outfile:
            outfile.write(vcf)


if __name__=="__main__":
        main()
