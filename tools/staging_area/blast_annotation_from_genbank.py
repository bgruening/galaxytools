#!/usr/bin/env python

from Bio import SeqIO
import sys, os

def index_genbank_features(gb_record, feature_type, qualifier) :
    answer = []
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            start = feature.location.start.position
            end = feature.location.end.position
            answer.append((start, end, feature.qualifiers))
            """
            if qualifier in feature.qualifiers :
                
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
            """
    return answer

def __main__():
    genbank_file = sys.argv[2]
    input_file = sys.argv[1]
    output_file = open(sys.argv[3], 'w')

    index = index_genbank_features(SeqIO.read(open(genbank_file,"r"), "genbank"), 'CDS', 'note')

    ret = []
    for line in open(input_file):
        values = line.strip().split('\t')
        if not values: continue
        orf_name, start, end = values[0], int(values[1]), int(values[2])

        for elem in index:
            if start >= elem[0] and end <= elem[1]:
                temp_hash = {'start': start, 'end': end}
                temp_hash.update(elem[2])
                ret.append(temp_hash)
                locus = elem[2].get('locus_tag', '')
                output_file.write('%s\t%s\t%s\t%s\t%s\n' % (orf_name, start, end, ' '.join(elem[2]['product']), ' '.join(locus) ))

if __name__ == "__main__" :
    __main__()


