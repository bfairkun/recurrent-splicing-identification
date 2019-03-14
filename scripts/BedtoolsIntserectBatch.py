#!/usr/bin/env python

import os
import sys
# USAGE:
# python BedtoolsIntersectBatch.py <FileIn>
#
# FileIn should be a tab-seperated file of <FileA> <FileB> <Output>. For each
# line in FileIn, this script will make a system call to a command constructed
# as follows:

# bedtools intersect -wa -a <FileA> -b <FileB> | bedtools intersect -wa -a - -b <FileC> > <Output>


FileIn = sys.argv[1]
if FileIn == "-":
    fh=sys.stdin
else:
    fh = open(FileIn, 'rU')

for line in fh:
    FileA, FileB, FileC, Output = line.strip('\n').split('\t')
    cmd = "bedtools intersect -wa -a {FileA} -b {FileB} | bedtools intersect -wa -a - -b {FileC} > {Output}".format(FileA=FileA, FileB=FileB, FileC=FileC, Output=Output)
    exit_code = os.system(cmd)
    sys.stderr.write(cmd + "\nfinished with exit code " + str(exit_code) + "\n")
fh.close()
