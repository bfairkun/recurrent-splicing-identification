# -*- coding: utf-8 -*-
#!/usr/bin/env python2

import os
import time
import re
import sys
import itertools
import cStringIO
import argparse
import distutils.dir_util
import gzip

def flush_stringio_buffers_to_disk(StringIOBufferDict, FilesToAppendInsteadOfWrite=None):
    """
    Given a dictionary that maps filepaths (the keys) to stringio objects (the values), this function will write out the stringio objects to the files, making directories as necessary if they are included in the filepath
    """
    for filename, stringio_buffer in StringIOBufferDict.items():
        distutils.dir_util.mkpath(os.path.dirname(filename))
        if filename in FilesToAppendInsteadOfWrite:
            fh = open(filename, 'a')
        else:
            fh = open(filename, 'w')
        stringio_buffer.seek(0)
        fh.write(stringio_buffer.read())
        fh.close()

def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()


def main():
    parser = argparse.ArgumentParser(usage = '''
    USAGE: python {myfile} -I [SnapTronDatabase.gz] -S [SampleKey.txt] [Params]
    '''.format(myfile=sys.argv[0])
    )
    parser.add_argument('-I','--InputFile', metavar = "[filepath]", help="REQUIRED. The input .gz snaptron file. Use 'stdin' or '-' to read from stdin", required=True)
    parser.add_argument('-S', '--SampleKey', metavar = "[filepath]", help="REQUIRED. two-column tab-delimited file that describes the snaptron sample_ids (usually a 5-ish digit code) and the corresponding path to write LeafCutter input files to", required=True)
    parser.add_argument("--quiet", help="OPTIONAL. quiet the output verbosity", action="store_true")
    parser.add_argument("--MaxLinesInMemory", metavar = "[integer]", help="OPTIONAL. Max number of output lines in memory files (buffer) at once before flushing to files on disk. Default = 100000000", type=int, default=100000000)
    parser.add_argument("--CreateEmptyOutputFiles", help="OPTIONAL. Create a file for every sample in the SampleKey file, regardless of if it ends up being an empty file without any corresponding junctions in the Snaptron file", action="store_true")
    args = parser.parse_args()

    # Make a dictionary of {snaptron_sampleid:outfilepath}
    snaptron_target_samples_dict = {}
    FilesOpened = set()
    with open(args.SampleKey) as f:
        for line in f:
            snaptron_sampleid,  outfilepath = line.strip('\n').split('\t')
            if args.CreateEmptyOutputFiles:
                    distutils.dir_util.mkpath(os.path.dirname(outfilepath))
                    open(outfilepath, 'w').close()
                    FilesOpened.add(outfilepath)
                    snaptron_target_samples_dict[snaptron_sampleid] = outfilepath

    MemoryFileBufferDict = {}
    LinesInMemoryFiles = 0
    BarcodeMatchedReadCount = 0
    BarcodeUnmatchedReadCount = 0 
    tic = time.clock()

    # if args.InputFile == 'stdin' or args.InputFile == '-':
    #     SnaptronFile = sys.stdin
    # else:
    #     SnaptronFile = args.InputFile

    SnaptronFile = args.InputFile
    with gzip.open(SnaptronFile) as f:
        for i,line in enumerate(f):
            if i:
                # Read in Snaptron database according to field names provided
                # on snaptron website
                junction_id, chrom, start, end, length, strand, annotated,leftmotif, rightmotif, left_annotated, right_annotated, samples_counts, total_count, cov_sum, cov_avg, cov_med, sourcedata_id=line.strip('\n').split('\t')

                # samples_count field takes the format of ,{sample}:{junction
                # counts},{sample2}:{junction counts}...
                for entry in samples_counts.strip(',').split(','):
                    sample_id, junction_count = entry.split(':')

                    if sample_id in snaptron_target_samples_dict:
                        # write out line to memory files stored in
                        # MemoryFileBufferDict values. Need to initialize
                        # dictionary entry if doesn't already exist

                        fileout = snaptron_target_samples_dict[sample_id]
                        if fileout not in MemoryFileBufferDict:
                            MemoryFileBufferDict[fileout] = cStringIO.StringIO()
                        MemoryFileBufferDict[fileout].write("{}\t{}\t{}\t.\t{}\t{}\n".format(chrom,start,end,junction_count, strand))
                        LinesInMemoryFiles += 1

                        if LinesInMemoryFiles >= args.MaxLinesInMemory:
                            # write out the memory files to disk and clear the
                            # memory files
                            flush_stringio_buffers_to_disk(MemoryFileBufferDict,FilesToAppendInsteadOfWrite=FilesOpened)
                            FilesOpened.update(MemoryFileBufferDict.keys())
                            MemoryFileBufferDict = {}
                            LinesInMemoryFiles = 0

            if not args.quiet:
                if i% 100000 == 0:
                        print(time.strftime('%X') + ' | ' + format((i)/float(1000000), '0.1f') + 'millions junctions processed')
                        restart_line()
        flush_stringio_buffers_to_disk(MemoryFileBufferDict,FilesToAppendInsteadOfWrite=FilesOpened)
    toc = time.clock()
    if not args.quiet:
        print('DONE')
        print('Elapsed time: ' + str(toc-tic))
        print('{i} junctions processed'.format(i=i))
main()

