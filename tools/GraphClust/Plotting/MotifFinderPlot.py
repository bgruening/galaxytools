from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import glob
import pandas as pd
import itertools
import seaborn as sns


def plot_bar(ranges, colors, orig_names, cluster_nums):
    fig, ax = plt.subplots()
    for i, k in enumerate(sorted(ranges.keys())):
        ax.broken_barh(ranges[k], (i-0.25, 0.5), facecolors=colors[k])

    ax.set_xlim(0)
    ax.set_xlabel('position in sequence')
    ax.set_yticklabels(['']+[k+'-'+orig_names[k] for k in sorted(ranges.keys())])
    ax.grid(True)
    fig.suptitle('Structure motif prediction\nRegions with same color are prediticted to have similar structures')
    # Add the legend
    patches = [mpatches.Patch(color=cluster_nums[lab], label=lab) for lab in sorted(cluster_nums)]
    ax.legend(handles=patches, loc='best')  # , bbox_to_anchor=(1, 0.5), loc='center left')
    plt.savefig("motif_plot.png")


def parse_clusters():
    currentdir_files = sorted(list(glob.glob('*')))
    print "currentdir_files are: ", currentdir_files
    cluster_files = sorted(list(glob.glob('RESULTS/*.cluster.all')))
    if len(cluster_files) == 0:
        raise RuntimeError('Expected cluster.all search path is empty:{}'.format(cluster_files))
    palette = itertools.cycle(sns.color_palette("Set2", len(cluster_files)))
    ranges = defaultdict(list)
    colors = defaultdict(list)
    orig_names = defaultdict(list)
    cluster_nums = defaultdict(list)
    for cluster_file in cluster_files:
        cluster_color = palette.next()
        df_cluster = pd.read_csv(cluster_file, sep='\s+', header=None)
        for irow, row in df_cluster.iterrows():
            seq, start, end, strand = row[0].split("#")
            ranges[seq].append((int(start), int(end)-int(start)+1))
            colors[seq].append(cluster_color)
            assert row[1] == 'RESULT'
            cluster_nums['cluster-{}'.format(row[2])] = cluster_color
            assert row[9] == 'ORIGHEAD'
            orig_names[seq] = row[10]
    return ranges, colors, orig_names, cluster_nums

my_ranges, my_colors, my_orig_names, my_cluster_nums = parse_clusters()
plot_bar(my_ranges, my_colors, my_orig_names, my_cluster_nums)
