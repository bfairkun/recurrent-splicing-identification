#!/usr/bin/env python3
import gzip
import sys
import pandas as pd
InGzFile, OutGzFile = sys.argv[1:]

df = pd.read_csv(InGzFile, index_col=0,  header=0, sep=' ')
print(df)
df.sort_index(axis=1, inplace=True)
df.to_csv(OutGzFile, sep=' ', compression='gzip')
