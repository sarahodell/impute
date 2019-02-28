#!/usr/bin/env python

import sys
import argparse
from subprocess import Popen,PIPE
import string

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description:
                                     Takes a vcf file and converts it to a csv format with founder as rows and marker genotypes as columns
                                     with nucleotide information encoded as A for reference allele and B for alternate alleles. Founder lines will be given letter codes based on the order that they are listed in the vcf file ("A" for the first sample and so on).
                                     This csv file is formatted for use with R/qtl2. It requires that bcftools be installed.
                                     """)
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args

def get_founders(vcf):
    """Gets a list of the samples in the vcf file"""
    process=Popen(["bcftools", "query","-l",vcf],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    print(stderr)
    founders=stdout.split('\n')[:-1]    
    return founders

def get_tmp(vcf):
    process=Popen(['bcftools','query','-f','%ID[\tGT=%GT]\n',vcf],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    info = stdout.split('\n')[:-1]
    print(stderr)
    return info

    
def get_foundergenofile():
    args=parse_args()
    txt='ind'
    geno={}
    founders=get_founders(args.infile)
    numf=len(founders)
    print(founders)
    f=''
    lcode=string.ascii_uppercase
    for s in range(1,numf+1):
        geno[s]=[founders[s-1]]
        f+='{0},{1}\n'.format(founders[s-1],lcode[s-1])
    print "Writing founder codes to FounderCodes.csv"
    with open('FounderCodes.csv','w') as ffile:
        ffile.write(f)
    info=get_tmp(args.infile)
    for i in info:
        split=i.split('\t')
        marker = split[0]
        txt+=','+marker
        count=1
        if 2 in split[1:]:
            for n in split[1:]:
                if './.' in n:
                    n='NA'
                elif '1' in n:
                    n='A'
                else:
                    n='B'
                geno[count].append(n)
                count+=1
        else:
            for n in split[1:]:
                if './.' in n:
                    n='NA'
                elif '0' in n:
                    n='A'
                else:
                    n='B'
                geno[count].append(n)
                count+=1
    for j in geno.keys():
        txt+='\n'+ ','.join(geno[j])
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)




if __name__ == "__main__":
    get_foundergenofile()


