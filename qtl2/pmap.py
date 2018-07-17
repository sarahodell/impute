#!/usr/bin/env python

"""
Generates a physical map file in csv format.
This csv file is formatted for use with R/qtl2
"""

import argparse
from subprocess import Popen,PIPE

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description: Generates a physical map file in csv format.                                                             
This csv file is formatted for use with R/qtl2""")
    parser.add_argument("infile",type=str,help="""The input vcf file containing the markers and physical positions""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args


def call_bcftools(vcf):
    process=Popen(["bcftools","query","-f","%ID,%CHROM,%POS\n",vcf],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    print(stderr)
    return stdout


def get_pmap():
    args=parse_args()
    txt='marker,chr,pos\n'
    stdout=call_bcftools(args.infile)
    markers=stdout.split('\n')[:-1]
    for line in markers:
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


