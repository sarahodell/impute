#!/usr/bin/env python

import sys
import argparse

def parse_args():
    """ -h for info on arguments
    """
<<<<<<< HEAD
    parser = argparse.ArgumentParser(description="""Program description:
                                     Takes a vcf file and converts it to a csv format with founder as columns and marker genotypes as rows
                                     with nucleotide information in IUPAC format.
                                     This csv file is formatted for use with R/qtl2
                                     """)
=======
    parser = argparse.ArgumentParser(description="""Program description: Takes a vcf file and converts it to a csv format with founder rows and marker genotypes as columns. Genotypes are assumed homozygous, with reference allele as A and alternate alleles as B.                                                             
This csv file is formatted for use with R/qtl2""")
>>>>>>> 461ea4af76df207712b0fb3aabdea80831872a58
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args


def get_foundergenofile():
    args=parse_args()
    txt='ind'
    samples={0:'A',1:'B',2:'C',3:'D',4:'E',5:'F',6:'G',7:'H',8:'I',9:'J',10:'K',11:'L',12:'M',13:'N',14:'O',14:'P'}
    index=0
    with open(args.infile,'r') as infile:
        tmp=''
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
    for s in samples:
        txt+='\n'+samples[s]
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)



if __name__ == "__main__":
    get_foundergenofile()


