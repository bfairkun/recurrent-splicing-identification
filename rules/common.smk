import pandas as pd
import os

snaptron_samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
# print(expand("expand this {filepaths}", filepaths = snaptron_samples.loc[:,"local_snaptron_file"]))
# print(snaptron_samples.at)
# df.loc[df['B'] == 3, 'A']

# Get cell based on other column value
# print("{mystr}".format( mystr = snaptron_samples.loc[snaptron_samples["local_snaptron_file"] == "/project2/yangili1/snaptron/SRA2_junctions_uncompressed.gz", "snaptron_wget_link"][0]))

# print(snaptron_samples.index)
# print(snaptron_samples.loc["SRA2", "local_snaptron_file"])
wildcard_constraints:
    snaptron_database_file = "|".join(snaptron_samples["local_snaptron_file"]),
    snaptron_metadata_file = "|".join(snaptron_samples["local_metadata_file"])
