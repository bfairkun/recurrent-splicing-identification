# Snakemake workflow: recurrent-splicing-identification

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/recurrent-splicing-identification.svg?branch=master)](https://travis-ci.org/snakemake-workflows/recurrent-splicing-identification)


## Authors

* Benjamin Fair (@bfairkun)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/recurrent-splicing-identification).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.


Install dependencies with conda.
```
conda env create --file envs/recurrent-splicing.yaml
conda activate leafcutter
```

Note that even though I named this environment `leafcutter`, it does not contain leafcutter software... I could not install leafcutter with conda. So you have to install leafcutter manually (see leafcutter documentation page) for R, and also have the leafcutter scripts accesible in PATH.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. Configure cluster settings in `cluster-config.json`

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"

or by executing the included sbatch script to execute the snakemake process from a cluster

    sbatch snakemake.sbatch

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
