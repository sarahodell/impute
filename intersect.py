#!/usr/bin/env python
"""
Calculates the interest of parent assignments between Simulated RIL and FILLIN Projection Alignment
"""

import sys
import argparse
from subprocess import PIPE,Popen,STDOUT
import pandas as pd
import os.path
import numpy as np

def arg_parse():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("actual",type=str,help="""Bed file of actual parental assignments""")
    parser.add_argument("pred",type=str,help="""Bed file of predicted parental assignments""")
    parser.add_argument("--verbose",type=str,help="""Print out full output, with per parent percentages""")
    args = parser.parse_args()
    return args


def get_total(pred):
    total=0
    with open(pred,'r') as full:
        for line in full:
            l = line.split('\t')
            start = int(l[1])
            end = int(l[2])
            total+= end-start
    return total


def percent_heterozygous(infile):
    df = pd.read_table(infile)
    fulltotal=sum(df['end']-df['start'])
    hetdf = df[df['donor1'] != df['donor2']]
    perc_het = sum(hetdf['end']-hetdf['start'])
    return round(float(perc_het)/fulltotal,3),fulltotal


def split_samples(infile,p=False,r=False):
    df = pd.read_table(infile)
    df = df[df['donor1']==df['donor2']]
    if p==True:
        seg='B'
        outtag='_pred.bed'
    else:
        seg='A'
        outtag='.bed'
    samples = df['sample'].unique()
    chroms = df['chr'].unique()
    parents=np.unique(df['donor1'])
    for i in samples:
        for c in chroms:
            count = 1
            txt=''
            locs = df[df['sample']==i]
            for index,row in locs.iterrows():
                start = row['start']
                end = row['end']
                txt+='{0}\t{1}\t{2}\t{3}\n'.format(c,start,end,row['donor1']+'_'+str(count)+seg)
                count+=1
        with open('tmp/'+i+'_'+str(c)+outtag,'w') as outfile:
            outfile.write(txt)
    if r==True:
        return samples,chroms,parents

    
def split_parents(p,infile,outfile):
    txt=''
    hasparent=False
    with open(infile,'r') as full:
        for line in full:
            if p in line.split('\t')[3]:
                txt+=line
                hasparent=True
    if hasparent:
        with open(outfile,'w') as pfile:
            pfile.write(txt)
    return hasparent

        
def call_bedtools(afile,bfile):
    process = Popen(['bedtools','intersect','-wao','-a',afile,'-b',bfile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr = process.communicate()
    print(stderr)
    return stdout


def get_overlap(stdout):
    overlap=[]
    lines = stdout.split('\n')
    for l in lines[:-1]:
        t = l.split('\t')
        o = int(t[-1])
        overlap.append(o)
    return sum(overlap)


def format_out(all_samples,samples,chroms):
    txt="Line\tTotal % Correct\tChrom\t% per Parent\n"
    for r in samples:
        total=all_samples[r]['total']
        tmp="{0}\t{1}\t".format(r,total)
        for c in chroms:
            tmp+='{0}'.format(c)
            parents = all_samples[r][c]["parent"]
            perc = all_samples[r][c]["perc_correct"]
            for i in range(len(parents)):
                tmp+='\t{0}:{1}'.format(parents[i],perc[i])
            txt+=tmp+'\n'
    txt+='\nTotal Percentage Correct: {0}\n'.format(all_samples['fulltotal'])
    txt+='Percentage of Homozygous Calls Correct: {0}\n'.format(all_samples['homozygoustotal'])
    txt+='Percentage of Heterozygous Calls: {0}\n'.format(all_samples['het'])
    return txt


def intersect(sample,c,parents):
    r_actual = 'tmp/{0}_{1}.bed'.format(sample,c)
    r_pred='tmp/{0}_{1}_pred.bed'.format(sample,c)
    total=get_total(r_pred)
    per_parent={"parent":[],"perc_correct":[],"right":0,"total":0}
    per_parent['total']=total
    for p in parents:
        afile='tmp/{0}_{1}_{2}.bed'.format(sample,c,p)
        bfile='tmp/{0}_{1}_{2}_pred.bed'.format(sample,c,p)
        if split_parents(p,r_actual,afile) and split_parents(p,r_pred,bfile):
            p_total = get_total(bfile)
            out = call_bedtools(afile,bfile)
            right = get_overlap(out)
            per_parent['parent'].append(p)
            per_parent['right']+=right
            per_parent['perc_correct'].append(round(float(right)/p_total,3))
        else:
            print('Parent {0} not shared between actual and predicted in sample {1} chr {2}\n').format(p,sample,c)
    return per_parent

def main():
    args = arg_parse()
    actual = args.actual
    pred = args.pred
    samples,chroms,parents=split_samples(actual,r=True)
    split_samples(pred,p=True)
    all_samples= {}
    pright=0
    ptotal=0
    for r in samples:
        all_samples[r]={}
        r_total=0
        r_right=0
        for c in chroms:
            all_samples[r][c] = intersect(r,c,parents)
            pright+=(all_samples[r][c]['right'])
            ptotal+=(all_samples[r][c]['total'])
            r_right+=(all_samples[r][c]['right'])
            r_total+=(all_samples[r][c]['total'])
        all_samples[r]['total']=round(float(r_right)/r_total,3)
    all_samples['het'],fulltotal=percent_heterozygous(pred)
    all_samples['fulltotal'] = round(float(pright)/fulltotal,3)
    all_samples['homozygoustotal']=round(float(pright)/ptotal,3)
    output = format_out(all_samples,samples,chroms)
    with open('intersect_output.txt','w') as outfile:
        outfile.write(output)


if __name__ == "__main__":
    main()
