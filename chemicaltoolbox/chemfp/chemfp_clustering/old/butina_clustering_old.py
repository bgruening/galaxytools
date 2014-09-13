#!/usr/bin/env python
"""
    Modified version of code examples from the chemfp project.
    http://code.google.com/p/chem-fingerprints/
    Thanks to Andrew Dalke of Andrew Dalke Scientific!
"""

import chemfp
import sys
import os
import tempfile

temp_file = tempfile.NamedTemporaryFile()
temp_link = "%s.%s" % (temp_file.name, 'fps')
temp_file.close()
os.system('ln -s %s %s' % (os.path.realpath(sys.argv[1]), temp_link) )


chemfp_fingerprint_file = temp_link
tanimoto_threshold = float(sys.argv[2])
outfile = sys.argv[3]
processors = int(sys.argv[4])


def get_hit_indicies(hits):
    return [id for (id, score) in hits]

out = open(outfile, 'w')
dataset = chemfp.load_fingerprints( chemfp_fingerprint_file )

chemfp.set_num_threads( processors )
search = dataset.threshold_tanimoto_search_arena(dataset, threshold = tanimoto_threshold)
#search = chemfp.search.threshold_tanimoto_search_symmetric (dataset, threshold = tanimoto_threshold)

# Reorder so the centroid with the most hits comes first.
# (That's why I do a reverse search.)
# Ignore the arbitrariness of breaking ties by fingerprint index
results = sorted( (  (len(hits), i, hits) for (i, hits) in enumerate(search.iter_indices_and_scores())  ),reverse=True)


# Determine the true/false singletons and the clusters
true_singletons = []
false_singletons = []
clusters = []

seen = set()

for (size, fp_idx, hits) in results:
    if fp_idx in seen:
        # Can't use a centroid which is already assigned
        continue
    seen.add(fp_idx)
    print size, fp_idx, hits
    if size == 1:
        # The only fingerprint in the exclusion sphere is itself
        true_singletons.append(fp_idx)
        continue

    members = get_hit_indicies(hits)
    # Figure out which ones haven't yet been assigned
    unassigned = [target_idx for target_idx in members if target_idx not in seen]

    if not unassigned:
        false_singletons.append(fp_idx)
        continue

    # this is a new cluster
    clusters.append( (fp_idx, unassigned) )
    seen.update(unassigned)

len_cluster = len(clusters)
#out.write( "#%s true singletons: %s\n" % ( len(true_singletons), " ".join(sorted(dataset.ids[idx] for idx in true_singletons)) ) )
#out.write( "#%s false singletons: %s\n" % ( len(false_singletons), " ".join(sorted(dataset.ids[idx] for idx in false_singletons)) ) )

out.write( "#%s true singletons\n" % len(true_singletons) )
out.write( "#%s false singletons\n" % len(false_singletons) )
out.write( "#clusters: %s\n" % len_cluster )

# Sort so the cluster with the most compounds comes first,
# then by alphabetically smallest id
def cluster_sort_key(cluster):
    centroid_idx, members = cluster
    return -len(members), dataset.ids[centroid_idx]

clusters.sort(key=cluster_sort_key)


for centroid_idx, members in clusters:
    centroid_name = dataset.ids[centroid_idx]
    out.write("%s\t%s\t%s\n" % (centroid_name, len(members), " ".join(sorted(dataset.ids[idx] for idx in members))))
    #ToDo: len(members) need to be some biggest top 90% or something ...

for idx in sorted(true_singletons):
    out.write("%s\t%s\n" % (dataset.ids[idx], 0))

out.close()
os.remove( temp_link )
