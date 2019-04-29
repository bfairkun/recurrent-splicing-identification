#!/usr/bin/env python3
import gzip
import sys
InGzFile, OutGzFile = sys.argv[1:]

with gzip.open(InGzFile, mode='rt') as f:
    with gzip.open(OutGzFile, 'wt') as fout:
        chromosome, *samples = f.readline().strip('\n').split()
        fout.write(chromosome + ' ' + ' '.join(sorted(samples)) + '\n')
        for line in f:
            junction, *counts = line.strip('\n').split()
            fout.write(junction + ' ' + ' '.join([x for _, x in sorted(zip(samples,counts))]) + '\n')
