#!/usr/bin/env python

import argparse
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

def arg_parse():
    parser=argparse.ArgumentParser(description="""Program description:
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.
1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::
    #chrom  start  end  donor1    donor2
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg
2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.

This code is modified from Ryan Dale (https://gist.github.com/daler/c98fc410282d7570efc3)

""")
    parser.add_argument("ideofile",type=str,help="""Tab-delimited file of breakpoint parental assignments""")
    parser.add_argument("out",type=str,default="image.png",help="""Name out image file to be outputed""")
    parser.add_argument("--colors",type=str,default=None,help="""tab-delimited File with parent and corresponding rbs color code, if None, colors are randomly generated""")
    parser.add_argument("--title",type=str,default="Parental Donors",help="Title of plot")
    parser.add_argument("--physical",type=bool,default=True,help="In units of Mb or cM? If true, than physical distance (Mb)")
    args = parser.parse_args()
    return args


def chromosome_collections(df, y_positions,height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width=False
    if 'width' not in df.columns:
        del_width=True
        df['width']=df['end'] - df['start']
    for s,group in df.groupby('sample'):
        print s
        yrange=(y_positions[s],height)
        xranges = group[['start','width']].values
        yield BrokenBarHCollection(
            xranges,yrange,facecolors=group['colors1'],**kwargs)
        yrange=(y_positions[s]-0.4,height)
        yield BrokenBarHCollection(
            xranges,yrange,facecolors=group['colors2'],**kwargs)
    if del_width:
        del df['width']


def make_rgb(parents):
    """Randomly generates rgb colors for each NAM parent for plotting 
    ideograms. Output is a dictionary with parent names as keys and 
    rgb lists as values"""
    color_lookup = {}
    for i in parents:
        rgb = numpy.random.random_sample(3)
        rgb = [round(j,2) for j in rgb]
        color_lookup[i]=rgb
    return color_lookup


def make_legend(parents,color_lookup):
    """Creates a legend in the plot with colors for each parent """
    legend = []
    for p in parents:
        pcolor=color_lookup[p]
        legend+=[Line2D([0],[0],color=pcolor,lw=4,label=p)]
    #legend+=[Line2D([0],[0],color=(0.0,0.0,0.0),lw=4,label='Cent')]
    return legend

def get_colors(colors,parents):
    color_lookup={}
    if colors==None:
        color_lookup=make_rgb(parents)
    else:
        with open(colors,'r') as infile:
            for line in infile:
                l = line[:-1].split('\t')
                color_lookup[l[0]]=l[1]
    return color_lookup


def draw_plot():
    args=arg_parse()
    
    # Read in file from output of make_ideogram()    
    ideo = pd.read_table(args.ideofile,skiprows=1,names=['sample','chr','start','end','donor1','donor2'])
    samples=list(ideo['sample'].unique())
    parents=ideo['donor1'].append(ideo['donor2']).unique()
    color_lookup=get_colors(args.colors,parents)
    # Height of each ideogram
    chrom_height=0.4
    # Spacing between consecutive ideograms
    chrom_spacing=1
    #Height of the gene track. Should be smaller than 'chrom_spacing' in order to fit correctly
    gene_height=0.2
    # Padding between the top of a gene track and its corresponding ideogram
    gene_padding=0.1
    # Width, height (in inches)
    figsize=(10,10)
    #Decide which chromosomes to use
    chromosome_list=samples

    # Keep track of the y positions for ideograms and genes for each chromsome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase={}
    gene_ybase={}
    chrom_centers={}
    
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase+= chrom_height + chrom_spacing

    #Filter out chromosomes not in our list
    #ideo = ideo[ideo.sample.apply(lambda x:x in chromosome_list)]
    ideo['end'] = pd.to_numeric(ideo['end'])
    ideo['start'] = pd.to_numeric(ideo['start'])
    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    #Add two new columns for the colors (same if homozygous, different if heterozygous)
    ideo['colors1'] = ideo['donor1'].apply(lambda x: color_lookup[x])
    ideo['colors2']=ideo['donor2'].apply(lambda x: color_lookup[x])
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    #Now all we have to do is call out function for the ideogram data
    for collection in chromosome_collections(ideo,chrom_ybase,chrom_height):
        ax.add_collection(collection)

    #Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')
    if args.physical==True:
        ax.set_xlabel('Position (Mb)')
    else:
        ax.set_xlabel('Position (cM)')
    ax.set_ylabel('Line Number')
    ax.set_title(args.title)


    box = ax.get_position()
    ax.set_position([box.x0,box.y0, box.width * 0.8, box.height])

    #put a legend to the right of the current axis
    legend_elements=make_legend(parents,color_lookup)
    ax.legend(handles=legend_elements,loc='center left',bbox_to_anchor=(1,0.5))

    plt.savefig(args.out)


if __name__=="__main__":
   draw_plot()               




