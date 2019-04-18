def local_snaptron_file_TO_snaptron_wget_link(wildcards):
    """
    Based on the database file (stored as a snakemake wildcard), return the corresponding wget link to download snaptron data, as found in the samples.tsv file (stored as snaptron_samples data frame)
    """
    return(snaptron_samples.loc[snaptron_samples["local_snaptron_file"] == wildcards.snaptron_database_file, "snaptron_wget_link"][0])

def local_metadata_file_TO_metadata_wget_link(wildcards):
    return(snaptron_samples.loc[snaptron_samples["local_metadata_file"] == wildcards.snaptron_metadata_file, "metadata_wget_link"][0])


rule download_snaptron_database:
    output:
        "{snaptron_database_file}"
    params:
        local_snaptron_file_TO_snaptron_wget_link
    shell:
        """
        wget -O {output} {params}
        """

rule download_snaptron_metadata:
    output:
        "{snaptron_metadata_file}"
    params:
        local_metadata_file_TO_metadata_wget_link
    shell:
        """
        wget -O {output} {params}
        """

rule download_chrom_sizes:
    output:
        chromesizes = "MiscData/hg38.chrome.sizes",
        chromosomal_genome_bed = "MiscData/bedfiles/ChromosomalGenome.bed"
    shell:
        """
        wget -q -O - https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes | awk -F '\\t' -v OFS='\\t' 'NR<=24' | tee {output.chromesizes} | awk -F '\\t' -v OFS='\\t' '{{ print $1, "1", $2 }}' > {output.chromosomal_genome_bed}
        """

rule make_junc_filelist:
    input:
        lambda wildcards: snaptron_samples.loc[{wildcards.sample}, "local_metadata_file"]
    output:
        "junction_filelist/{sample}.tsv"
    shell:
        """
        awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $1, "juncfiles/{wildcards.sample}/"$1".junc"}}' {input} > {output}
        """

rule make_junc_files:
    """
    This rule will make all the junction files for each sample in
    junction_filelist/{sample}.tsv. Since snaptron contains tens of thousands
    of samples, rather than add all of these as target files to the DAG (which
    might make snakemake painfully slow to build the DAG), just specify a log
    file as output that gets modified after all junction files are made.
    """
    input:
        junction_filelist = "junction_filelist/{sample}.tsv",
        snaptron_file = lambda wildcards: snaptron_samples.loc[{wildcards.sample}, "local_snaptron_file"]
    output:
        "logs/make_junc_files/{sample}.log"
    shell:
        """
        python scripts/SplitSnaptronBySampleId.py -I {input.snaptron_file} -S {input.junction_filelist} --CreateEmptyOutputFiles > {output} && echo "job done" >> {output}
        """

rule make_bed_file_for_each_chrom:
    input:
        "MiscData/hg38.chrome.sizes"
    output:
        "MiscData/bedfiles/{chrom}.bed"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$1=="{wildcards.chrom}" {{print $1, 1, $2}}' {input} > {output}
        """

