#!/usr/bin/env python

"""
Generates a physical map file in csv format.
This csv file is formatted for use with R/qtl2
"""

import sys
import argparse

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description: Generates a physical map file in csv format.                                                             
This csv file is formatted for use with R/qtl2""")
    parser.add_argument("infile",type=str,help="""The input intermediate csv file generated by bcftools query""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args

def get_pmap():
    args=parse_args()
    txt='marker,chr,pos\n'
    with open(args.infile,'r') as infile:
        for line in infile:
            info = line.split(',')
            marker = info[0]
            chrom=info[1]
            pos=float(info[2])/1e6
            tmp='{0},{1},{2}\n'.format(marker,chrom,pos)
            txt+=tmp
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)



if __name__ == "__main__":
    get_pmap()


