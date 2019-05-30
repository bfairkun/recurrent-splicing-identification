# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
include: "rules/common.smk"

rule all:
    input:
        "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.numers.gz".format(Samples_TargetJunctions=Samples_TargetJunctions),
         # "leafcutter/differential_splicing/{leafcutter_outprefix}_effect_sizes.txt".format(leafcutter_outprefix = leafcutter_ds_outprefix),
        # "leafcutter/differential_splicing/{leafcutter_outprefix}_cluster_significance.txt".format(leafcutter_outprefix = leafcutter_ds_outprefix),
        "logs/delete_temp_files/{Samples_TargetJunctions}.log".format(Samples_TargetJunctions=Samples_TargetJunctions),
        "logs/GatherDyanamicBed/{Samples_TargetJunctions}.log".format(Samples_TargetJunctions=Samples_TargetJunctions),
        # expand("{snaptron_databses}", snaptron_databses = snaptron_samples["local_snaptron_file"]),
        # expand("MiscData/TargetBedFilesPerChrom/{chrom}.bed", chrom=["chr%d"%x for x in range(1,23)]+['chrX','chrY']),
        # expand("MiscData/make_BedtoolsIntersectBatch_input_files/{chrom}", chrom=["chr%d"%x for x in range(1,23)]+['chrX','chrY']),
        # expand("logs/make_juncfile_per_chrom/{chrom}", chrom=["chr%d"%x for x in range(1,23)]+['chrX','chrY']),
        # expand("leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/leafcutter_perind.counts.gz", chrom=["chr%d"%x for x in range(1,23)]+['chrX','chrY']),

        # expand("{snaptron_databses}", snaptron_databses = snaptron_samples["local_metadata_file"]),
        # expand("junction_filelist/{sample}.tsv", sample=snaptron_samples.index),
        # expand ( "logs/make_junc_files/{sample}.log", sample=snaptron_samples.index),
        # "junction_filelist/AllSamples",
        # "leafcutter_out/clustering/headtest/leafcutter_perind.counts.gz"
        # "leafcutter_out/clustering/AllSamples/leafcutter_perind.counts.gz"

include: "rules/GetDataForLeafcutter.smk"
include: "rules/leafcutter.smk"
include: "rules/other.smk"
