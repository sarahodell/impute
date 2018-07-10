#!/usr/bin/env python


import argparse
import pandas as pd
import numpy as np
from scipy import stats

def arg_parse():
    parser=argparse.ArgumentParser(description="""Program description: Uses the Ogut 2015 maize genetic map to approximate the genetic distance between two physical positions""")
    parser.add_argument("chrom",type=int,help="""Chromosome number""")
    parser.add_argument("start",type=float,help="""Start position (bp)""")
    parser.add_argument("end",type=float,help="""End position (bp)""")
    parser.add_argument("--ref",type=str,default="v4",help="""Reference version of physical positions (either v2, v3, or v4) """)
    args = parser.parse_args()
    return args


def approx_cM(chrom,start,end,ref='v4'):
    """Arguments:
    chr: (int) chromosome number 1..10
    start: (int) physical bp position of start
    end: (int) physical bp position of end (end>start)
    ref (optional): (str) reference version of physical positions (either v2,v3, or v4)

    Returns:
    The approximate genetic distance between the start and end positions based on the Ogut 2015 genetic map (int)
    """
    ogutpath='/group/jrigrp/Share/annotations/genetic_map/ogut2015/ogut_fifthcM_map_agp{0}.txt'.format(ref)
    ogutmap=pd.read_table(ogutpath,sep='\s+',header=None,
                          names=['SNP_ID','SNP_newID','chr','pos','cM'])
    ogutmap.dropna(axis=0,inplace=True)
    ogutmap['chr']=pd.to_numeric(ogutmap['chr'])
    ogutmap['pos']=pd.to_numeric(ogutmap['pos'])
    #Grab the chromosome
    cmap=ogutmap[ogutmap['chr']==chrom]
    cmap.reset_index(inplace=True)
    #Closest ogut map markers above and below the start and end postions
    below=cmap.iloc[cmap[cmap['pos']<=start]['pos'].idxmax(),]
    above=cmap.iloc[cmap[cmap['pos']>=end]['pos'].idxmin(),]
    xpoints=[below['pos'],above['pos']]
    ypoints=[below['cM'],above['cM']]
    #Linear regression of genetic position on physical position using two closest points
    slope,intercept,r_value,p_value,stderr = stats.linregress(xpoints,ypoints)
    start_cM=slope*start + intercept
    end_cM=slope*end + intercept
    #return the genetic distance
    return abs(end_cM-start_cM)



if __name__ == "__main__":
    args=arg_parse()
    distance=approx_cM(args.chrom,args.start,args.end,ref=args.ref)
    print("The approximate genetic distance between {0} and {1} is {2} cM ({3} cM/bp)".format(int(args.start),int(args.end),round(distance,3),round(distance/(args.end-args.start),6)))
