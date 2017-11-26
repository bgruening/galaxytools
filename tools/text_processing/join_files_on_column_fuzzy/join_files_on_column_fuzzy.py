#!/usr/bin/env python

import os
import argparse
import sys

def main(args):

    if args.header:
        h1 = True
        h2 = True
    else:
        h1 = False
        h2 = False

    cache = list()
    out = open(args.outfile, 'w+')
    write_buffer = list()

    def _readline(header = False):
        with open(args.f2) as handle2:
            for line in handle2:
                line = line.strip()
                if header:
                    header = False
                    yield line
                    continue
                if not line:
                    continue
                columns = line.split(args.sep)
                value2 = columns[args.c2-1]
                yield columns, float(value2)

    def fill_cache():
        try:
            cache.append(next(it))
        except StopIteration:
           pass

    it = _readline(header = h2)

    with open(args.f1) as handle1:
        for line in handle1:
            line = line.strip()
            if h1:
                h1 = False
                seconda_header = next(it)
                if args.add_distance:
                    out.write('%s\t%s\t%s\n' % (line, seconda_header, args.unit))
                else:
                    out.write('%s\t%s\n' % (line, seconda_header))
                continue
            if not line:
                continue
            columns = line.split(args.sep)
            value1 = float(columns[args.c1-1])
            _cache = list()
            fill_cache()
            while cache:
                _c, value2 = cache.pop(0)
                upper_bound = value1 + args.distance
                if args.unit == 'absolute':
                    if value2 <= upper_bound and value2 >= (value1 - args.distance):
                        line_template = '%s\n'
                        abs_dist = abs(value1 - value2)
                        if args.add_distance:
                            line_template = '%s\t' + str(abs_dist) + '\n'
                        write_buffer.append([abs_dist, line_template % '\t'.join( columns + _c )])
                        _cache.append([_c, value2])
                        fill_cache()
                    elif value2 > upper_bound:
                        # if the value from list 2 is bigger then the current value, he will be taken into the next round
                        _cache.append([_c, value2])
                    elif value2 < upper_bound:
                        # if the value from list 2 is smaller then the currecnt value, check the next one of list 2
                        fill_cache()
                elif args.unit == 'ppm':
                    ppm_dist = abs((value1 - value2) / value1 * 1000000)
                    if ppm_dist <= args.distance:
                        line_template = '%s\n'
                        if args.add_distance:
                            line_template = '%s\t' + str(ppm_dist) + '\n'
                        write_buffer.append([ppm_dist, line_template % '\t'.join( columns + _c )])
                        _cache.append([_c, value2])
                        fill_cache()
                    elif ppm_dist > args.distance:
                        _cache.append([_c, value2])
                    elif ppm_dist < args.distance:
                        fill_cache()
            if args.closest and write_buffer:
                write_buffer.sort(key=lambda x: x[0])
                out.write(write_buffer[0][1])
            else:
                for _dist, line in write_buffer:
                    out.write(line)
            write_buffer = list()
            cache = _cache
    out.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Merge two files on a common column the fuzzy way.')
    parser.add_argument('--f1', required=True)
    parser.add_argument('--f2', required=True)
    parser.add_argument('--c1', type=int, required=True, help="Column in file 1 to be merged on.")
    parser.add_argument('--c2', type=int, required=True, help="Column in file 2 to be merged on.")
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--header', action='store_true', help="The files have a header line at the beginning.")
    parser.add_argument('--closest', action='store_true', help="Only report the closest match.")
    parser.add_argument('--add_distance', action='store_true', help="Add addional column with the distance between the two values.")
    parser.add_argument('--sep', type=str, default="\t", help="Files are separated by this separator.")
    parser.add_argument('--distance', type=float, default="0.2", help="Maximal allowed distance.")
    parser.add_argument('--unit', choices=['ppm', 'absolute'], default='absolute')
    args = parser.parse_args()

    main(args)


