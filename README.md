Snakemake workflow: rna-count-salmon

This is the [Snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) workflow for RNA-Seq read count, powered by [Salmon](https://salmon.readthedocs.io/en/latest/). Optional quality controls are performed by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregated by [MultiQC](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)

Each tool belong to their respective authors.

# Installation

## Install conda

In order to install conda, please follow the guide lines described on the official [Conda web-page](https://docs.continuum.io/anaconda/install/).

This will likely be:

```{sh}
# Conda pre-requisites
apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

# Download latest conda installer
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh

# Run the installer
bash ~/Downloads/Anaconda3-2019.03-Linux-x86_64.sh
```

You can alternatively try the script: `setup.sh` within the `.circleci ` directory, by typing `bash .circleci/setup.sh`. But remember, this is provided for machine without either conda, nor miniconda installed.

## Install Snakemake pre-requisites within Conda

We are going to use Python 3.7. No other options. Forget about python 2.7, it will be dropped out of support before 2020'. In [rna-count-salmon](https://github.com/gustaveroussy/rna-count-salmon), we use [formatted strings](https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting) available at python 3.6. It is not excluded that further development will include python 3.7 exclusive behaviors. Go for latest python !

A [yaml]() description of all dependencies is embedded with this work flow under `envs/workflow.yaml`. Install all dependencies with:

```{sh}
conda env create --file envs/workflow.yaml
```

In which `python`, `datrie`, `snakemake`, `pandas` and `git` are non-optional and essential. `Pytest`, is used to perform unit-testing. Finally, `jinja2`, `pygraphviz`, `flask`, `openssl`, `zlib`, and `networkx` are here to handle dynamic snakemake's reporting.

Every single time you will use this pipeline, you have to activate that conda environment. The activation command will not be repeated everywhere in this wiki, so please, don't forget it!

## Clone this workflow

The final step to use this workflow is to clone it to your computer. This can be easily done with git (which you have installed in your virtual environment few seconds ago!

Here we go:

```{sh}
git clone https://bitbucket.org/tdayris/rna-count-salmon.git
```

# Testing

## Unit testing

Scripts embedded within this pipeline can be tested with `pytest`. Go within the tests sub-directory and run:

```{sh}
make all-unit-tests
```

# Scripts sanlity testing

Separate tests can be run on each script with:

```{sh}
make config-tests
make design-tests
make aggregation-tests
```

## Smoke testing

The sanity tests are used for continuous integration, just go within the tests sub-directory and run:

```{sh}
make conda-tests
make ci-tests
```

# Configuration

## The general config file (yaml)

This is a [yaml](https://en.wikipedia.org/wiki/YAML) file. It could also be [json](https://en.wikipedia.org/wiki/JSON) formatted, however its name must not be changed. Yaml format has been preferred since
users requested a more human readable configuration file to write.

Note that the scripts `prepare_config.py` creates a custom configuration file from command line. Try it!

This configuration file contains the following sections (in any orders):

### ref

This is a very simple section: two values with evident significations. We need two paths, one to the fasta formatted genome sequence, and one other to the gtf formatted genome annotation. By GTF formatted, I mean GTF. No GFF, no GFF3, no GFF2, no GTF.gz, no GTF.bz2, nothing else than GTF. The same goes with the fasta file.

Example:
```{yaml}
ref:
  fasta: path/to/fasta.fa
  gtf: /path/to/gtf.gtf
```

These paths can be either relative or absolute.

### params

This section contains additional command line arguments. One for all, **do not try to modify threading options** here. Do not try to modify  neither logging, nor temporary directories. These options are al
ready handled.

We have the following values:

copy_extra: Extra parameters for bash cp
salmon_index_extra: Extra parameters for kallisto index rule
salmon_quant_extra: Extra parameters for kallisto quant rule

Example:
```{yaml}
params:
  copy_extra: --parents --verbose
  libType: ISF
  salmon_index_extra: \-\-kmerLen 5
  salmon_quant_extra: \-\-noBiasLengthThreshold \-\-minAssignedFrags 1 \-\-noEffectiveLengthCorrection
    \-\-noLengthCorrection \-\-fasterMapping \-\-noFragLengthDist \-\-allowDovetail
    \-\-numPreAuxModelSamples 0 \-\-numAuxModelSamples 0
```

### workflow

This section is used to activate or deactivate sections of the pipeline. In fact, if you just want to quantify and no quality control (you may have done them aside of the [rna-count-salmon](https://github.com/gustaveroussy/rna-count-salmon) pipeline), then turn these flag to `false` instead of `true`.

Warning, it is case sensitive!

Example:
```
workflow:
  fastqc: true
  multiqc: true
  aggregate: true
```

### General

The following parameters do not belong to any section, let them be at the top level of the yaml file:

design: The path to the design file
workdir: The path to the working directory
threads: The maximum number of threads used
singularity_docker_image: The image used within singularity
cold_storage: A list of cold storage mount points

Example:
```
design: design.tsv
workdir: .
threads: 1
singularity_docker_image: docker://continuumio/miniconda3:4.4.10
cold_storage:
  - /media
```

### Conclusion

A complete config.yaml file would look like this:

```
design: design.tsv
workdir: .
threads: 1
singularity_docker_image: docker://continuumio/miniconda3:4.4.10
cold_storage:
  - /media
ref:
  fasta: /path/to/genome/sequence.fa
  gtf: /path/to/genome/annotaton.gtf
workflow:
  fastqc: true
  multiqc: true
  aggregate: true
params:
  copy_extra: --parents --verbose
  libType: ISF
  salmon_index_extra: \-\-kmerLen 5
  salmon_quant_extra: \-\-noBiasLengthThreshold \-\-minAssignedFrags 1 \-\-noEffectiveLengthCorrection
    \-\-noLengthCorrection \-\-fasterMapping \-\-noFragLengthDist \-\-allowDovetail
    \-\-numPreAuxModelSamples 0 \-\-numAuxModelSamples 0
```

Remember, there is a python script here to build this file from command line!

## The design file (tsv)

At this stage, we need a [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) file describing our analysis.

Note that the script `prepare_design.py` tries to creates this file for you. Due to the very wide range of sample naming conventions, this script may not produce expected results: it relies on alphabetical order to gather *fastq pairs*.

It must contain the following columns:

* Sample_id: the name of each samples
* Upstream_file: path to the upstream fastq file

The optional columns are:

* Downstream_file: path to downstream fastq files (usually your R2 in paired-end libraries)
* Any other information

Remember, there is a python script here to build this file from command line!

# Usage

## Local usage

Here is your command line:

```{sh}
snakemake
```

Yes. Nothing more as long as you have [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/), and [Salmon](https://salmon.readthedocs.io/en/latest/) in you `PATH`.

However, conda provides an efficient way to dynamically build virtual environments and to double check your tool versions, that's why we recommend:

```{sh}
snakemake --use-conda
```

And if you want to be sure that your OS is not interfering with your run, then use:

```{sh}
snakemake --use-conda --singularity
```

I won't detail more of the Snakemake possibilities, please see [their documentation](https://snakemake.readthedocs.io/en/stable/executable.html). However, many users like the following:

```{sh}
snakemake --use-conda --singularity --reason --printshellcmds
```

With these supplementary arguments, you will have the executed command lines printed to your STDOUT, and each rule will give you the reason why it is being executed. Nice, isn't it ?

## PBS Torque

In addition to the previous arguments, one might want to use this pipeline on a PBS Torque cluster.

```{sh}
snakemake --cluster " qsub -l walltime=00:{resources.time_min}:00,nodes=1:ppn={threads},mem={resources.mem_mb}mb -V " --use-conda --singularity --reason --printshellcmds --jobs 20  --restart-times 3
```

## Slurm

In addition to the previous arguments, one might want to use this pipeline on a Slurm cluster.

```{sh}
snakemake--cluster " sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time={resources.time_min} " --use-conda --singularity --reason --printshellcmds --jobs 20  --restart-times 3
```


# Pipeline content

## Global workflow

This workflows takes fastq files, genome sequences and annotations as input, and returns abundance estimates along side with optional quality metrics.

If you use this pipeline, cite them all, please!

## FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) stands for FastQ Quality Control. There is no publication associated with this tool, however, it remains an inescapable classic in bioinformatics.

This tool compiles a lot of informations about the raw reads obtained from sequencers. Actually, this tool does not perform any process **required** for the whole splicing analysis; however, it's a good practice.

Citation:

* Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010).

## MultiQC

[MultiQC](https://multiqc.info/), just like FastQC, do not have any other purpose than quality metrics. It gathers all Flagstat and all FastQC individual metrics into one single report.

Citation:

* Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.


## Salmon

[Salmon](https://salmon.readthedocs.io/) is a tool for transcript quantification from RNA-seq data. It uses pseudo-mapping to compute quantification estimates on transcripts.

Citation:

* Patro, Rob, et al. “Salmon provides fast and bias-aware quantification of transcript expression.” Nature Methods (2017).

## Snakemake

[Snakemake](https://snakemake.readthedocs.io) is a pipeline/workflow manager written in python. It is used to handle the tools interaction, dependencies, command lines and cluster reservation. It is the skeleton of this pipeline. This pipeline is powered by the [Snakemake-Wrappers](https://snakemake-wrappers.readthedocs.io), the [Snakemake Workflows](https://github.com/snakemake-workflows), and the [conda](https://anaconda.org) project.

Citation:

* Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
