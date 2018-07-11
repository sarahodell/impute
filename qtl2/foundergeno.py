#!/usr/bin/env python

import sys
import argparse

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description:
                                     Takes a vcf file and converts it to a csv format with founder as rows and marker genotypes as columns
                                     with nucleotide information encoded as A for reference allele and B for alternate alleles.
                                     This csv file is formatted for use with R/qtl2
                                     """)
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    parser.add_argument("founders",type=str,help="""File with a list of founders in the order they are listed in the vcf file""")
    args=parser.parse_args()
    return args


def get_foundergenofile():
    args=parse_args()
    txt='ind'
    samples={}
    key=0
    with open(args.founders,'r') as ffile:
        for line in ffile:
            samples[key]=line[:-1]
            key+=1
    with open(args.infile,'r') as infile:
        for line in infile:
            info = line.split(',')
            marker = info[0]
            txt+=',{0}'.format(marker)
            count=0
            for i in info[3:]:
                if './.' in i:
                    n='NA'
                elif '0' in i:
                    n='A'
                else:
                    n='B'
                samples[count]+=',{0}'.format(n)    
                count+=1
    for k in range(key):
        txt+='\n'+samples[k]
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)



if __name__ == "__main__":
    get_foundergenofile()


