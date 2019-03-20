rule prepare_juncfiles_for_clustering:
    input:
        bedfile = junction_intersect_bed,
        chromosome_beds = "MiscData/bedfiles/{chromosome}.bed",
        make_junc_files = make_junc_files_log,
        juncfile = juncfiles
    output:
        leafcutter_juncfile = "leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfile.{{chromosome}}.txt".format(Samples_TargetJunctions=Samples_TargetJunctions),
        scratch_dir_touch = config["scratch_prefix"] + "leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/{{chromosome}}/touch".format(Samples_TargetJunctions=Samples_TargetJunctions),
    log:
        "logs/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/{{chromosome}}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        # mkdir -p leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/{wildcards.chromosome}/
        touch {output.scratch_dir_touch}
        awk -F'\\t' -v OFS='\\t' '{{  print $1,  "{input.bedfile}", "{input.chromosome_beds}" ,"{config[scratch_prefix]}leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/{wildcards.chromosome}/" gensub("/",".","g",$1)  }}' {juncfiles} | python scripts/BedtoolsIntserectBatch.py - 2> {log}
        awk -F'\\t' -v OFS='\\t' '{{ print   "{config[scratch_prefix]}leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfiles/{wildcards.chromosome}/" gensub("/",".","g",$1) }}' {juncfiles} > {output.leafcutter_juncfile}
        """

rule leafcutter_cluster:
    input:
        leafcutter_juncfile = "leafcutter/prepare_juncfiles_for_clustering/{Samples_TargetJunctions}/juncfile.{{chromosome}}.txt".format(Samples_TargetJunctions=Samples_TargetJunctions),
    output:
        config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.gz".format(Samples_TargetJunctions=Samples_TargetJunctions),
        config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind_numers.counts.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    log:
        "logs/leafcutter_cluster/{Samples_TargetJunctions}/{{chromosome}}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    params:
        rundir = config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        leafcutter_cluster.py -j {input.leafcutter_juncfile} -r {params.rundir} &> {log}
        """

rule rearrange_leafcutter_cluster_counts:
    input:
        config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    output:
        config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.rearranged.gz".format(Samples_TargetJunctions=Samples_TargetJunctions),
    shell:
        """
        python3 scripts/RearrangeColumnsAndSplitFractionsOfLeafcutterClusterCountTable.py {input} {output}
        """

rule make_numers_and_denoms:
    input:
        config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.rearranged.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    output:
        numers = config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.numers".format(Samples_TargetJunctions=Samples_TargetJunctions),
        denoms = config["scratch_prefix"] + "leafcutter/clustering/{Samples_TargetJunctions}/{{chromosome}}/leafcutter_perind.counts.denoms".format(Samples_TargetJunctions=Samples_TargetJunctions),
    shell:
        """
        zcat {input} | perl -lne 'if ($.==1) {{print}} else {{$_ =~ s/\d+\///g; print}}' > {output.denoms}
        zcat {input} | perl -lne 'if ($.==1) {{print}} else {{$_ =~ s/\/\d+//g; print}}' > {output.numers}
        """

rule merge_numers_and_denoms:
    input:
        numers_to_gather = numers_to_gather,
        denoms_to_gather = denoms_to_gather
    output:
        # merged count tables. Include a header, for easier reading into R
        numers_merged = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.numers.gz".format(Samples_TargetJunctions = Samples_TargetJunctions),
        denoms_merged = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.denoms.gz".format(Samples_TargetJunctions = Samples_TargetJunctions),

        # leafcutter differential splicing analysis scripts expect a header without an entry for the row names.
        numers_merged_for_leafcutter_analysis = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind_numers.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        awk 'NR==1 {{print}} FNR>1 {{print}}' {input.numers_to_gather} | gzip - > {output.numers_merged}
        awk 'NR==1 {{print}} FNR>1 {{print}}' {input.denoms_to_gather} | gzip - > {output.denoms_merged}
        awk 'NR==1 {{printf $2; for (i=3; i <= NF; i++) printf FS$i; print NL}} FNR>1 {{print}}' {input.numers_to_gather} | gzip - > {output.numers_merged_for_leafcutter_analysis}
        """

rule delete_temp_files:
    """
    Used a rule to do this rather than the temp() snakemake directive for
    temporary files because leafcutter produces a file for each junction file.
    So using a rule with 'rm -rf' avoids adding all those potentially thousands
    of files to the DAG
    """
    input:
        numers_merged_for_leafcutter_analysis = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind_numers.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    output:
        "logs/delete_temp_files/{Samples_TargetJunctions}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        "rm -rf {temporary_clusterfiles} 2> {output}"
