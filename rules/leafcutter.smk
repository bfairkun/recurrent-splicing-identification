# If juncfile for leafcutter_cluster is given in config, use it, and don't
# check that the .junc count files exist (To avoid creating DAGS that may have
# tens of thousands of files)
if config["leafcutter_cluster_juncfiles"]:
    leafcutter_cluster_juncfiles = config["leafcutter_cluster_juncfiles"]
    make_junc_files_log = []
else:
    leafcutter_cluster_juncfiles = "junction_filelist/AllSamples"
    make_junc_files_log = expand ( "logs/make_junc_files/{sample}.log", sample=snaptron_samples.index),

rule leafcutter_cluster:
    input:
        juncfile = leafcutter_cluster_juncfiles,
        make_junc_files_log = make_junc_files_log
        # Rather than explicitly requiring that .junc files listed in juncfile
        # actually exist, just require the log file for all make_junc_files
        # output as input. This keeps the DAG from getting cluttered with tens
        # of thousands of files.
    output:
        "leafcutter_out/clustering/{junctionsfile}/leafcutter_perind.counts.gz".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
        "leafcutter_out/clustering/{junctionsfile}/leafcutter_perind_numers.counts.gz".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
    log:
        "logs/leafcutter_cluster/{junctionsfile}.log".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
    params:
        rundir = "leafcutter_out/clustering/{junctionsfile}/".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) )
    shell:
        """
        leafcutter_cluster.py -j {input.juncfile} -r {params.rundir} &> {log}
        """

rule leafcutter_cluster_perchrom:
    input:
        juncfile = "MiscData/make_leafcutter_cluster_juncfiles/{chrom}",
        make_junc_files_log = "logs/make_juncfile_per_chrom/{chrom}"
        # Rather than explicitly requiring that .junc files listed in juncfile
        # actually exist, just require the log file for all make_junc_files
        # output as input. This keeps the DAG from getting cluttered with tens
        # of thousands of files.
    output:
        "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/leafcutter_perind.counts.gz",
        "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/leafcutter_perind_numers.counts.gz",
    log:
        "logs/leafcutter_cluster_perchrom/{chrom}.log",
    params:
        rundir = "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/"
    shell:
        """
        leafcutter_cluster.py -j {input.juncfile} -r {params.rundir} &> {log}
        """
# rule Make_perind_numers_file:

# rule Make_perind_denom_file:

