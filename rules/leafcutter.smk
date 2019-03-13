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

if config["leafcutter_cluster_by_chrom"]:
    ToGather = expand("leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfile.{chromosome}.txt", Samples_TargetJunctions=Samples_TargetJunctions, chromosome=Chromosome_list)
else:
    ToGather = expand("leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfile.{chromosome}.txt", Samples_TargetJunctions=Samples_TargetJunctions, chromosome="ChromosomalGenome")



rule prepare_juncfiles_for_clustering:
    input:
        bedfile = junction_intersect_bed,
        chromosome_beds = "MiscData/bedfiles/{chromosome}.bed",
        make_junc_files = make_junc_files_log,
        juncfile = juncfiles
    output:
        leafcutter_juncfile = "leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfile.{{chromosome}}.txt".format(Samples_TargetJunctions=Samples_TargetJunctions),
    log:
        "logs/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/{{chromosome}}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        mkdir -p leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles
        awk -F'\\t' -v OFS='\\t' '{{ print $1,  "{input.bedfile}", "{input.chromosome_beds}" ,"leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/" $1  }}' {juncfiles} | python scripts/BedtoolsIntserectBatch.py - 2> {log}
        awk -F'\\t' -v OFS='\\t' '{{ "leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/" $1  }}' {juncfiles} > {output.leafcutter_juncfile}
        """

rule gather_prepares:
    input:
        leafcutter_juncfile = ToGather
    output:
        "Gathered"

# rule leafcutter_cluster:
#     input:
#         juncfile = leafcutter_cluster_juncfiles,
#         make_junc_files_log = make_junc_files_log
#         # Rather than explicitly requiring that .junc files listed in juncfile
#         # actually exist, just require the log file for all make_junc_files
#         # output as input. This keeps the DAG from getting cluttered with tens
#         # of thousands of files.
#     output:
#         "leafcutter_out/clustering/{junctionsfile}/leafcutter_perind.counts.gz".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
#         "leafcutter_out/clustering/{junctionsfile}/leafcutter_perind_numers.counts.gz".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
#     log:
#         "logs/leafcutter_cluster/{junctionsfile}.log".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) ),
#     params:
#         rundir = "leafcutter_out/clustering/{junctionsfile}/".format(junctionsfile = os.path.basename(leafcutter_cluster_juncfiles) )
#     shell:
#         """
#         leafcutter_cluster.py -j {input.juncfile} -r {params.rundir} &> {log}
#         """

# rule leafcutter_cluster_perchrom:
#     input:
#         juncfile = "MiscData/make_leafcutter_cluster_juncfiles/{chrom}",
#         make_junc_files_log = "logs/make_juncfile_per_chrom/{chrom}"
#         # Rather than explicitly requiring that .junc files listed in juncfile
#         # actually exist, just require the log file for all make_junc_files
#         # output as input. This keeps the DAG from getting cluttered with tens
#         # of thousands of files.
#     output:
#         "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/leafcutter_perind.counts.gz",
#         "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/leafcutter_perind_numers.counts.gz",
#     log:
#         "logs/leafcutter_cluster_perchrom/{chrom}.log",
#     params:
#         rundir = "leafcutter_out/clustering/AllSamples_PerChrom/{chrom}/"
#     shell:
#         """
#         leafcutter_cluster.py -j {input.juncfile} -r {params.rundir} &> {log}
#         """
# # rule Make_perind_numers_file:

# # rule Make_perind_denom_file:

