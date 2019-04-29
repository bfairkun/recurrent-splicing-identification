#!/usr/bin/env python2

import argparse
import sys
import gzip
import os
from collections import defaultdict
from operator import itemgetter

parser = argparse.ArgumentParser(description="Make aggregated bed file of junction counts from leafcutter numer or denom gzipped count table. Useful for visualizing differential splicing raw data in IGV.")
parser.add_argument('-C','--CountsFileIn',  help="leafcutter counts.gz table. Be careful, header should contain entry for row names", required=True)
parser.add_argument('-G', '--GroupsFileIn', help="leafcutter groups file", required=True)
parser.add_argument('-O', '--OutputFilePrefix', default="", help="Output file prefix")
parser.add_argument('-OD', '--OutputDirectory', default="./", help="Output directory")
parser.add_argument("--quiet", help="OPTIONAL. quiet the output verbosity", action="store_true")
args = parser.parse_args()

if not os.path.exists(args.OutputDirectory):
        os.makedirs(args.OutputDirectory)

Groups = defaultdict(list)
with open(args.GroupsFileIn, 'rU') as f:
    for line in f:
        sample, group = line.strip('\n').split('\t')
        Groups[group].append(sample)

with gzip.open(args.CountsFileIn, 'rU') as f:
    header = f.readline().strip('\n').split()

for group in Groups:
    if not args.quiet: print("Processing {}...".format(group))
    group_indices = [i for i, val in enumerate(header) if val in set(Groups[group])]
    with open(args.OutputDirectory + args.OutputFilePrefix + group + ".bed", "w") as fout:
        with gzip.open(args.CountsFileIn, 'rU') as f:
            for i,line in enumerate(f):
                if i>=1:
                    l = line.strip('\n').split()
                    chrom, start, stop, cluster = l[0].split(':')
                    group_count_sum = (sum([int(i) for i in itemgetter(*group_indices)(l)]))
                    fout.write("{}\t{}\t{}\t{}\t{}\t.\n".format(chrom, start, stop, l[0], group_count_sum))
