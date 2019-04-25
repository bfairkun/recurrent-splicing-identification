rule prepare_juncfiles_for_clustering:
    """Intersects junction files at junction_intersect_bed defined in config
    file, so that we don't have to run leafcutter_cluster on whole genome if it
    is not necessary. This alalso handles splitting the junction files by
    chromosome to run leafcutter chromosome by chromosome for better
    parallelization on cluster and also so that memory requirements are more
    reasonable."""
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
    log:
        "logs/rearrange_leafcutter_cluster_counts/{Samples_TargetJunctions}/{{chromosome}}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        python3 scripts/ReorderLeafcutterClusterColumns.py {input} {output} 2> {log}
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
        numers_merged_for_leafcutter_analysis = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.numers.gz".format(Samples_TargetJunctions=Samples_TargetJunctions)
    output:
        "logs/delete_temp_files/{Samples_TargetJunctions}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        "rm -rf {temporary_clusterfiles} 2> {output}"

rule AggregateGroupJunctionCounts:
    """Creates bed files for easy IGV visualization of the aggregate of raw
    junction counts of the two groups definied in the groups file (each group
    contains many samples)"""
    input:
        numers_merged = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.numers.gz".format(Samples_TargetJunctions = Samples_TargetJunctions),
        denoms_merged = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind.counts.denoms.gz".format(Samples_TargetJunctions = Samples_TargetJunctions),
        groupfile = config["leafcutter_groupfile"]
    output:
        numers_bed = dynamic("leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/numers.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions)),
        denoms_bed = dynamic("leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/denoms.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions)),
    shell:
        """
        python scripts/MakeJunctionBed.py -C {input.numers_merged} -G {input.groupfile} -O leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/numers.
        python scripts/MakeJunctionBed.py -C {input.denoms_merged} -G {input.groupfile} -O leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/denoms.
        """

rule AggregateGroupJunctionRatioCounts:
    input:
        numers_bed = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/numers.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions),
        denoms_bed = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/denoms.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions),
    output:
        ratios_bed = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/ratio.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions)
    shell:
        """
        paste {input.numers_bed} {input.denoms_bed} | awk -F'\\t' -v OFS='\\t' '$11>0 {{ print $1,$2,$3,$4,$5/$11,$6  }} $11==0 {{ print $1,$2,$3,$4,$5,$6 }}' > {output.ratios_bed}
        """

rule GatherDynamicBed:
    input:
        numers_bed = dynamic("leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/numers.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions)),
        denoms_bed = dynamic("leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/denoms.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions)),
        ratios_bed = dynamic("leafcutter/clustering/{Samples_TargetJunctions}/Merged/AggregateJuncFiles/ratio.{{group}}.bed".format(Samples_TargetJunctions = Samples_TargetJunctions))
    output:
        "logs/GatherDyanamicBed/{Samples_TargetJunctions}.log".format(Samples_TargetJunctions=Samples_TargetJunctions)
    shell:
        """
        echo 'done' > {output}
        """


rule leafcutter_ds:
    input:
        numers_merged_for_leafcutter_analysis = "leafcutter/clustering/{Samples_TargetJunctions}/Merged/leafcutter_perind_numers.gz".format(Samples_TargetJunctions=Samples_TargetJunctions),
        groupfile = config["leafcutter_groupfile"]
    output:
        "leafcutter/differential_splicing/{leafcutter_outprefix}_effect_sizes.txt".format(leafcutter_outprefix = leafcutter_ds_outprefix),
        "leafcutter/differential_splicing/{leafcutter_outprefix}_cluster_significance".format(leafcutter_outprefix = leafcutter_ds_outprefix),
    threads: 4
    params:
        outputprefix = "-o leafcutter/differential_splicing/{leafcutter_outprefix}".format(leafcutter_outprefix = leafcutter_ds_outprefix),
        exons = "-e R_project/gencode.v26.exons.txt.gz"
    log:
        "logs/leafcutter_ds/{leafcutter_outprefix}.log".format(leafcutter_outprefix = leafcutter_ds_outprefix)
    shell:
        """
        mkdir -p leafcutter/differential_splicing/
        leafcutter_ds.R -p {threads} {params.outputprefix} {params.exons} {input.numers_merged_for_leafcutter_analysis} {input.groupfile} &> {log}
        """

