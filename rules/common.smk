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

Chromosome_list = ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"]

wildcard_constraints:
    snaptron_database_file = "|".join(snaptron_samples["local_snaptron_file"]),
    snaptron_metadata_file = "|".join(snaptron_samples["local_metadata_file"])

# If juncfile for leafcutter_cluster is given in config, use it, and don't
# check that the .junc count files exist (To avoid creating DAGS that may have
# tens of thousands of files)
if config["juncfiles"]:
    juncfiles = config["juncfiles"]
    make_junc_files_log = []
else:
    juncfiles = "junction_filelist/AllSamples"
    make_junc_files_log = expand ( "logs/make_junc_files/{sample}.log", sample=snaptron_samples.index),

if config["junction_intersect_bed"]:
    junction_intersect_bed = config["junction_intersect_bed"]
else:
    junction_intersect_bed = "MiscData/bedfiles/ChromosomalGenome.bed"

Samples_TargetJunctions = os.path.basename(juncfiles) + "_" + os.path.basename(junction_intersect_bed)
temporary_clusterfiles = config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/".format(Samples_TargetJunctions=Samples_TargetJunctions)

if config["leafcutter_cluster_by_chrom"]:
    numers_to_gather = expand(config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{chromosome}/leafcutter_perind.counts.numers", Samples_TargetJunctions=Samples_TargetJunctions, chromosome=Chromosome_list)
    denoms_to_gather = expand(config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{chromosome}/leafcutter_perind.counts.denoms", Samples_TargetJunctions=Samples_TargetJunctions, chromosome=Chromosome_list)
else:
    numers_to_gather = expand(config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{chromosome}/leafcutter_perind.counts.numers", Samples_TargetJunctions=Samples_TargetJunctions, chromosome="ChromosomalGenome")
    denoms_to_gather = expand(config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{chromosome}/leafcutter_perind.counts.denoms", Samples_TargetJunctions=Samples_TargetJunctions, chromosome="ChromosomalGenome")
