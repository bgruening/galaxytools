#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import glob
import pandas as pd
import itertools
import seaborn as sns
import numpy as np
import sys
from matplotlib import cm

#from mpltools import style
#from mpltools import layout
#import matplotlib.pyplot as plt
# plt.style.use(['ggplot','seaborn-whitegrid'])
import re 

def sorted_alphanum( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
def plot_bar(ranges, colors, orig_names, cluster_nums):
    fig = plt.figure(figsize=(10,len(ranges)*0.6))
    ax = fig.add_subplot(211)
    sorted_ranges_keys = sorted(ranges.keys(),reverse=True, key=lambda s: int(s.replace('SEQ','')))

    for i, k in enumerate(sorted_ranges_keys):
        ranges_odd=[r for ri, r in enumerate(ranges[k]) if ri%2==1]
        ranges_even=[r for ri, r in enumerate(ranges[k]) if ri%2==0]

        colors_odd = [c for ci, c in enumerate(colors[k]) if ci%2==1]
        colors_even = [c for ci, c in enumerate(colors[k]) if ci%2==0]
        ax.broken_barh(ranges_odd, (i,0.5), facecolors=colors_odd,alpha=0.8)
        ax.broken_barh(ranges_even, (i-0.3,0.5), facecolors=colors_even,alpha=0.8)

    ax.set_xlim(0)
    ax.set_xlabel('position in sequence')
    ax.set_yticks(np.arange(-1, len(ranges)))
    hide_SEQ = True
    if hide_SEQ:
        ax.set_yticklabels(['']+[orig_names[k] for k in sorted_ranges_keys])
    else:
        ax.set_yticklabels(['']+[k+'-'+orig_names[k] for k in sorted_ranges_keys])
    ax.grid(True)
    # fig.suptitle('Structure motif prediction\nRegions with same color are prediticted to have similar structures')
    # Add the legend
    print(cluster_nums)
    patches = [mpatches.Patch(color=cluster_nums[lab], label=lab) for lab in sorted_alphanum(cluster_nums.keys())]
    # ax.legend(#handles=patches, loc='best',     bbox_to_anchor=(0, 0, 0.5, 0.5), bbox_transform=fig.transFigure)#bbox_to_anchor=(1.2, 1.05))#, loc='center left')
    import math

    ax.legend(handles=patches,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=8,ncol=math.ceil(len(patches)/30.0))

    fig.savefig("motif_plot.png", bbox_inches='tight')
    fig.savefig("motif_plot.svg", bbox_inches='tight')


def parse_clusters(res_dir):
    res_dir = glob.escape(res_dir)
    currentdir_files = sorted(list(glob.glob('*')))
    print ("currentdir_files are: ", currentdir_files)
    print ("RESULTS_files are: ", sorted(list(glob.glob(res_dir+'/*'))))

    cluster_files = sorted(list(glob.glob(res_dir+'/*.cluster.all')))
    if len(cluster_files) == 0:
        raise RuntimeError('Expected cluster.all search path is empty:{}'.format(cluster_files))
    colors_ggplot_colorblind = ["#BBBBBB",  "#56B4E9", "#E69F00","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    colors_Vega = cm.Vega20(np.linspace(0, 1.0, len(cluster_files)))

    use_colors = colors_ggplot_colorblind
    if len(cluster_files) > len(colors_ggplot_colorblind):
    	use_colors = colors_Vega
    palette = itertools.cycle(sns.color_palette(use_colors
        #"Paired", max(len(cluster_files),5)
        ))

    ranges = defaultdict(list)
    colors = defaultdict(list)
    orig_names = defaultdict(list)
    cluster_nums = defaultdict(list)
    for cluster_file in cluster_files:
        cluster_color = next(palette)
        df_cluster = pd.read_csv(cluster_file, sep='\s+', header=None)
        for irow, row in df_cluster.iterrows():
            seq, start, end, strand = row[0].split("#")
            ranges[seq].append((int(start), int(end)-int(start)+1))
            colors[seq].append(cluster_color)
            assert row[1] == 'RESULT'
            cluster_nums['cluster-{}'.format(row[2])] = cluster_color
            assert row[9] == 'ORIGHEAD'
            if(len(row)<11):
                print('WARNING: row {}, missing input label, using fasta id instead:\n{}\n'.format(irow, ' '.join(row)))
                orig_names[seq] = row[8]
            else:
                orig_names[seq] = row[10]
    return ranges, colors, orig_names, cluster_nums

if (len(sys.argv)>1):
    results_dir = sys.argv[1]
else:
    results_dir = "RESULTS/"
my_ranges, my_colors, my_orig_names, my_cluster_nums = parse_clusters(results_dir)
plot_bar(my_ranges, my_colors, my_orig_names, my_cluster_nums)
