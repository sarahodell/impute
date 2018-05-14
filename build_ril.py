#!/usr/bin/env python

from subprocess import Popen, PIPE
import pandas as pd
import sys
import argparse
import numpy as np

def get_args():
    parser=argparse.ArgumentParser(description="""Program description: Create mixed vcf file from multiple donor vcf files. Optional: make a certain proportion of data in the output vcf file missing, for imputation purposes.""")
    parser.add_argument("infile",type=str,help="""The input file: outpute file of write_rils() function in SimulateRILs """)
    parser.add_argument("outfile",type=str,help="""Name of output vcf file""")
    parser.add_argument("--drop",type=bool,help=""""Boolean: drop data? """)
    parser.add_argument("--droprate",type=float,help="""Float, proportion of data to drop(i.e. 0.2)""")
    args = parser.parse_args()
    if args.drop:
        print "Drop Frequency set to {0}".format(args.droprate)
    return args
        
def main():
    args=get_args()
    bedfile = pd.read_table('{0}'.format(args.infile),sep='\t')
    samples = bedfile.columns[3:]
    header = ''
    process = Popen(['bcftools','view','-h','c10_hmp31_edit_founders.vcf.gz'],stdout=PIPE,stderr=PIPE)
    stdout,stderr = process.communicate()
    header+=stdout
    for sample in samples:
        vcf = header
        for index,row in bedfile.iterrows():
            chrom = row[0]
            start = row[1]
            end = row[2]
            donorfile = '{0}_c10_hmp321.vcf.gz'.format(row[sample])
            query = '{0}:{1}-{2}'.format(chrom,start,end)
            process = Popen(['bcftools','view','-H','-r',query,donorfile],stdout=PIPE,stderr=PIPE)
            stdout,stderr = process.communicate()
            if args.drop == True:
                new_out = ''
                for line in stdout.split('\n'):
                    draw = np.random.random_sample()
                    if draw >= args.droprate:
                        new_out+=line
                        new_out+='\n'
                vcf+=new_out
            else:
                vcf+=stdout
        with open('{0}_{1}'.format(sample,args.outfile),'w') as outfile:
            outfile.write(vcf)

    
if __name__=="__main__":
        main()
