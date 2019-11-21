[![GitHub actions status](https://github.com/tdayris-perso/rna-count-salmon/workflows/CI/badge.svg?branch=master)](https://github.com/tdayris-perso/rna-count-salmon/workflows/CI/badge.svg?branch=master)

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

### FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) stands for FastQ Quality Control. There is no publication associated with this tool, however, it remains an inescapable classic in bioinformatics.

This tool compiles a lot of informations about the raw reads obtained from sequencers. Actually, this tool does not perform any process **required** for the whole splicing analysis; however, it's a good practice.

Citation:

* Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010).

### MultiQC

[MultiQC](https://multiqc.info/), just like FastQC, do not have any other purpose than quality metrics. It gathers all Flagstat and all FastQC individual metrics into one single report.

Citation:

* Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.


### Salmon

[Salmon](https://salmon.readthedocs.io/) is a tool for transcript quantification from RNA-seq data. It uses pseudo-mapping to compute quantification estimates on transcripts.

Citation:

* Patro, Rob, et al. “Salmon provides fast and bias-aware quantification of transcript expression.” Nature Methods (2017).

### Snakemake

[Snakemake](https://snakemake.readthedocs.io) is a pipeline/workflow manager written in python. It is used to handle the tools interaction, dependencies, command lines and cluster reservation. It is the skeleton of this pipeline. This pipeline is powered by the [Snakemake-Wrappers](https://snakemake-wrappers.readthedocs.io), the [Snakemake Workflows](https://github.com/snakemake-workflows), and the [conda](https://anaconda.org) project.

Citation:

* Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

## Understanding methodology

If you want to understand the whole ideas behind this pipeline, please read the following (tools above are not repeated):

1. Roberts, Adam, et al. [Improving RNA-Seq expression estimates by correcting for fragment bias.](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-3-r22) Genome Biology 12.3 (2011): 1.
1. Love, Michael I., Hogenesch, John B., Irizarry, Rafael A. [Modeling of RNA-seq fragment sequence bias reduces systematic errors in transcript abundance estimation.](https://www.nature.com/articles/nbt.3682) Nature Biotechnology 34.12 (2016).
1. Soneson C, Love MI and Robinson MD. [Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](file:///home/tdayris/Documents/Projects/B17_FROL/Treatment_modif/report.html) F1000Research 2016, 4:1521
1. Srivastava A, Sarkar H, Gupta N, Patro R; [RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985), Bioinformatics, Volume 32, Issue 12, 15 June 2016, Pages i192–i200
1. Bray N.L. et al. . ( 2016) [Near-optimal probabilistic RNA-seq quantification.](https://www.nature.com/articles/nbt.3519) Nature Biotech., 34(5), 525-527.
1. Ceppellini, r., Siniscalco, M. & Smith, C.A. [The estimation of gene frequencies in a random-mating population](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-1809.1955.tb01360.x) Ann. Hum. Genet. 20, 97–115 (1955)
1. Dempster, A.P., Laird, N.M. & rubin, D.B. J. R. [Maximum Likelihood from Incomplete Data via the EM Algorithm](https://www.jstor.org/stable/2984875) Stat. Soc. Ser. B 39, 1–38 (1977)
1. Chambers, John M., and Trevor J. Hastie, eds. [Statistical models in S.](https://sksznwbat01.storage.googleapis.com/MDQxMjA1MzAxMg==01.pdf) Vol. 251. Pacific Grove, CA: Wadsworth & Brooks/Cole Advanced Books & Software, 1992.
1. Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, [Differential analysis of RNA-Seq incorporating quantification uncertainty](https://www.nature.com/articles/nmeth.4324), Nature Methods (2017)
1. Kanitz, A., Gypas, F., Gruber, A. J., Gruber, A. R., Martin, G., & Zavolan, M. (2015). [Comparative assessment of methods for the computational inference of transcript isoform abundance from RNA-seq data.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0702-5) Genome biology, 16(1), 150.
1. Dillies, M. A., Rau, A., Aubert, J., Hennequet-Antier, C., Jeanmougin, M., Servant, N., … & Guernec, G. (2013).[ A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis.](https://academic.oup.com/bib/article/14/6/671/189645) Briefings in bioinformatics, 14(6), 671-683.
1. Storey, J. D., & Tibshirani, R. (2003). [Statistical significance for genomewide studies.](http://www.pnas.org/content/100/16/9440.short) Proceedings of the National Academy of Sciences, 100(16), 9440-9445.

# Frequently asked question by my fellow bioinformaticians on this pipeline

## What is cold storage?

On most clusters, a difference is being made between repositories designed to IO-intensive computations (hot storage), and repositories designed to archive (cold storage). Locally, system administrators allow (or not) the users to compute on archive-designed repositories.

Yet, copying all your fastq files at once before starting to work is time consuming and requires two copies of all your fastq files at the beginning of your work. Let Snakemake handle both copy creation and removal. This will save you time!

## I have seen a typo/coding error.

You should not. I am very happy to improve myself, this code, and to help you get more confidence in yourself. You have found a mistake? Open an issue. CHEERS and THANK YOU!

## I want to change/add/replace anything in your code. Can I?

You can open an issue, or fork this git repository, and make your own path. You are likely to have very good reasons to do so. Who knows, maybe you'll make a pull request?

# Frequently asked questions by my fellow biologists on this pipeline

## Understand FastQC content:

This file contains multiple graphs which will be described one after each other below. If you have any issue to understand a graph in spite of these details, please contact your bioinformatician.

### Basic Statistics

With basic statistics you can see a table containing the following:

* **Filename**: the name of the fastq file analyzed bu FastQC.
* **File type**: the type of file analyzed by FastQC. Usually, this will be "Conventional base calls".
* **Encoding**: the [encoding](https://en.wikipedia.org/wiki/FASTQ_format#Encoding) (a.k.a the method) used by the sequencer to write down quality metrics in the fastq file. Usually, this will be: "Sanger / Illumina 1.9".
* **Total Sequences**: the total number of reads sequenced. If this number is below your expectetions toward the sequencing platform, please contact both your bioinformatician and the sequencing platform.
* **Sequences flagged as poor quality**: Number of reads that were discarded. Usually, this will be 0 or very low.
* **Sequence length**: the length of the sequenced reads. If this number does not match with your expectations, please contact both your bioinformatician and the sequencing platform.
* **%GC**: the percentage of [Guanine (G)](https://en.wikipedia.org/wiki/Guanine) and [Cytosine (C)](https://en.wikipedia.org/wiki/Cytosine) in your sequenced reads. For humans, this is likely to be between 45% and 55%. For mice, this is usually a little below: 42% to 52%. If this value is far over, or far below this number, and your do not understand why, then contact both your bioinformatician and the sequencing platform.

### Per base sequence quality

This is a boxplot.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the quality of each base.

When a sequencer sequences a base, it attributes a quality on each measurement. This quality is called a [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score).

| Phred Quality Score | Probability of incorrect base call | Base call accuracy | Comment                    |
|:--------------------|:-----------------------------------|:-------------------|:---------------------------|
| 10                  | 1 in 10                            | 90%                | Not good                   |
| 20                  | 1 in 100                           | 99%                | Fair                       |
| 30                  | 1 in 1000                          | 99.9%              | Good                       |
| 40                  | 1 in 10000                         | 99.99%             | Very good                  |
| 50                  | 1 in 100000                        | 99.999%            | So good I've never seen it |

You can see in the background, three colored areas:

* **green**: with Phred scores over 28, the best ones
* **orange**: with Phred scores over 20, the fair ones
* **red**: with Phred scores below 20, the worst ones

On the graph, you should see:

* **Yellow boxplots**: The [quartiles](https://en.wikipedia.org/wiki/Quartile) of the distribution of the qualities over your reads at a given position.
* **Black moustaches**: The [deciles](https://en.wikipedia.org/wiki/Decile) of the distribution of the qualities over your reads at a given position.
* **Red line**: The mean quality of the distribution of the qualities over your reads at a given position.
* **Blue line**: The general mean of the quality of your reads

In RNA-Seq, quality scores usually are lower in two sections:

* *The first ten to fifteen bases on the read*: This is due to the (not-so) random priming used to keep your reads on your tile. Thermodynamically, all your reads will not hybridize as quick and as well as every other ones. This is what causes this lower Phred call.
* *The last bases of your read*: The longer the DNA polymerase keeps adding nucleotides, the higher is the probability that it will include wrong bases. Thus it also has higher probability to stop polymerizing. This leads to a lower Phred call. The longer your reads are, the stronger this effect will be. Usually, with 75bp long reads, your are likely to barely see this effect. With 250bp long reads, the effect will be clearly present.

For any other low quality values that you do not understand, please contact your bioinformatician.

### Per tile sequence quality

This is a heatmap.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the position on your tile. The heatmap color scheme goes from blue (good) to red (bad).

This pipeline is intended to quantify reads for RNA-Seq differential expression. Please, do not focus on this graph: errors within reads usually are not a problem for differential expression analysis. Greens and red strays over of your heatmap, random red dots, or colored circles/ellipses are usually not problematic.

If you are in doubt, please contact your bioinformatician.

### Per sequence quality scores

This is a classic line plot.

Horizontally (a.k.a the x axis), you have the mean Phred score. Vertically (a.k.a in y axis), you have the times this phred score was counted. This graph shows the count of bases called with a given Phred Score. Basically, we wan most of our counts around higher Phred scores.

See above, in [Per base sequence quality]() to have more information about Phred scores.

### Per base sequence content

This a four-overlaying line plots.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the percentage of usage of each base. The four lines refers to their respective bases (A, T, G, and C).

Usually, in RNA-Seq, you will see two/three strange things:

* *The first ten to fifteen bases on the read*: This is due to the (not-so) random priming used to keep your reads on your tile. Thermodynamically, all your reads will not hybridize as quick and as well as every other ones. You can see the sequence that was Thermodynamically selected during the sequencing process.
* *There is a higher proportion of As in the end of my reads*: Hooray! This is an artifact of Poly-A selection in RNA-Seq.
* *There is a higher proportion of Gs in the end of my reads*: Hooray! This is an artifact due to too short reads, leading the sequencer to add Gs. When it does not receive signal, sequencers usually add Gs.

Those three bias are perfectly fair in RNA-Seq and do not lead to errors when taken into account by your bioinformatician.

### Per sequence GC content

This is a two-overlaying line graph.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the percentage of usage of GCs. The blue line refers to a theoretical [(normal) distribution](https://en.wikipedia.org/wiki/Normal_distribution). The red line corresponds to the observed distribution.

If you are working on Humans/mice, you are likely to see the red and blue line overlying each other quite perfectly: means of both should be very close, as well as their standard-deviation. You should see only one continuous curve in either blue and red line.

If you see a two peaks, it usually signs a contamination.

### Per base N content

This is a classic line plot.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the percentage of Ns. N is called when sequencer could not call a base (A, T, C or G) with enough confidence. Basically, we do not want any N any where.

### Sequence Length Distribution

This is a classic line plot.

Horizontally (a.k.a the x axis), you have the length of your sequences. Vertically (a.k.a in y axis), you have the counts of reads having this given length. If this number does not match your expectations, please contact your bioinformatician and the sequencing platform.

### Sequence Duplication Levels

This is a two-overlaying line graph.

Horizontally (a.k.a the x axis), you have the count of this sequence duplication. Vertically (a.k.a in y axis), you have the count of the times a sequence has been seen this many time.

In RNA-Seq, we do want duplicated sequences, as the protocol is quantitative.

### Overrepresented sequences

This is a table giving the list of over-represented sequences and their proportion. Feel free to [blast](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) these sequences on NCBI.

The less overrepresented sequences, the better. However, this is usually not problematic in RNA-Seq, as long as it expected by the protocol.

### Adapter Content

This is a four-overlaying line plots.

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the percentage of usage of each adapter. The four lines refers to their respective databases (type of adapters).

We expect the as less adapters as possible.

In RNA-Seq, you are likely to see an increasing percentage of adapters in the end of your reads. This is usually not problematic and highlights shorter reads in your libraries. Contact your bioinformatician for more information.

### Kmer Content

You do know what a k-mer is. A dimer is a k-men of length 2. A trimer is a k-mer of length 3. A k-mer is a sequence of length *k*. Here, this graph shows heptamers (k-mers of length 7).

Horizontally (a.k.a the x axis), you have the position in the read in base pairs. Vertically (a.k.a in y axis), you have the percentage of usage of each kmer. Usually, this information is intended for bioinformaticians themselves. We barely use it in RNA-Seq, yet FastQC produces this graph.

You should not have to worry about it. It is useful in genome assembly, in meta-omics, etc.


## Understand MultiQC content:

### Step 1: Understand the content of FastQC

In this pipeline, most of the content of MultiQC is the one of FastQC. There are only one table and one line-plot added through Salmon. Have a read at the sections above!

### Fragment length distribution

Salmon, while quasi mapping (see below for more information on quasi mapping), has defined best positions of a read over the genome/transcriptome. If you have used paired-end sequencing protocol to sequence your RNA-Seq, then a pair of reads makes a "fragment".

The length of this fragment should match with your expectation according to the sequencing platform. If not, please contact your bioinformatician and the sequencing platform.

### Summary table

This summary table contains quite explicit information relative to the number of reads, number of pseudo-mapped reads, percentage of GC, and number of mapped reads.


## Why the `--seqbias` in Salmon ?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#seqbias), we can read that Salmon is able to learn and correct sequence-specific bias. It models the (not-so-)random hexamer printing bias and tries to correct it. By default, Salmon learns the sequence-specific bias parameters using 1,000,000 reads from the beginning of the input. Salmon uses a variable-length Markov Model (VLMM) to model the sequence specific biases at both the 5’ and 3’ end of sequenced fragments. This methodology generally follows that of [Roberts et al.](file:///home/tdayris/Documents/Projects/B17_FROL/Treatment_modif/report.html#SeqBiasArticle), though some details of the VLMM differ.

This options keeps us from trimming the ~13 first nucleotides in the beginning of every reads in our RNA-Seq data. Usually, these nucleotides present a lot of errors and/or enrichment (the preferential sequencing of fragments starting with certain nucleotide motifs).

## What is bootstrapping ?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#numbootstraps), we can read that Salmon is able to perform [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29). We strongly recommend to perform at least 35 bootstraps rounds in order to have a correct estimation of the technical variance within each of your samples.

Bootstrapped abundance estimates are done by resampling (with replacement) from the counts assigned to the fragment equivalence classes, and then re-running the optimization procedure, either the EM or VBEM (see below for more information), for each such sample. The values of these different bootstraps allows us to assess technical variance in the main abundance estimates we produce. Such estimates can be useful for downstream (e.g. differential expression) tools that can make use of such uncertainty estimates. The more samples computed, the better the estimates of varaiance, but the more computation (and time) required.

## Why the `--gcbias` in Salmon ?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#gcbias), we can read that Salmon is able to learn and correct fragment-level GC biases in the input data. Specifically, this model will attempt to correct for biases in how likely a sequence is to be observed based on its internal GC content.

A strong GC bias is often a sign of an enrichment in a given sequence, a contamination or a sequencing error. GC Bias correction does not impair quantification for samples without GC bias. For samples with moderate to high GC bias, correction for this bias at the fragment level has been shown to reduce isoform quantification errors: [Pato et al.](https://www.nature.com/articles/nmeth.4197), and [Love et al.](https://www.nature.com/articles/nbt.3682).

This bias is distinct from the primer biases previously described. Though these biases are distinct, they are not completely independent.

## How was the kmer length chosen in Salmon?

A k-mer is a word of length “k”. A Dimer is a k-mer of length 2. A Pentamer is a k-mer of length 5, etc. To be specific of a region, due to the length of the human genome, we strongly recommend a value above 27 (usually, the default value of 31 is enough, see [Srivastava et al.](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985)).

However, shorter reads are discarded. This means that any sequence of length below 31 will not appear here, or will present potential analysis bias.

## How was the sequencing library parameter chosen?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/library_type.html), we can read that all theoretical positions of reads are taken into account by Salmon.

There are numerous library preparation protocols for RNA-seq that result in sequencing reads with different characteristics. For example, reads can be single end (only one side of a fragment is recorded as a read) or paired-end (reads are generated from both ends of a fragment). Further, the sequencing reads themselves may be un-stranded or strand-specific. Finally, paired-end protocols will have a specified relative orientation.

* I = inward
* O = outward
* M = matching
* S = stranded
* U = un-stranded
* F = read 1 (or single-end read) comes from the forward strand
* R = read 1 (or single-end read) comes from the reverse strand
* A = Automatic detection

## What is mapping validation in Salmon ?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#validatemappings), we can read that Salmon can perform mapping validation:

One potential artifact that may arise from alignment-free mapping techniques is spurious mappings. These may either be reads that do not arise from some target being quantified, but nonetheless exhibit some match against them (e.g. contaminants) or, more commonly, mapping a read to a larger set of quantification targets than would be supported by an optimal or near-optimal alignment.

Salmon will run an extension alignment dynamic program on the quasi-mappings it produces. Moreover, Salmon makes use of an intelligent alignment cache to avoid re-computing alignment scores against redundant transcript sequences (e.g. when a read maps to the same exon in multiple different transcripts).

Quality trimming is more important before processing reads using this option.

## What is quasi-mapping?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-quasi-index-and-fmd-index-based-modes), we can read that:

One of the novel and innovative features of Salmon is its ability to accurately quantify transcripts using quasi-mappings. Quasi-mappings are mappings of reads to transcript positions that are computed without performing a base-to-base alignment of the read to the transcript. Quasi-mapping is typically much faster to compute than traditional (or full) alignments, and can sometimes provide superior accuracy by being more robust to errors in the read or genomic variation from the reference sequence.

The quasi-mapping was introduced in the work of [Srivastava et al](https://academic.oup.com/bioinformatics/article/32/12/i192/2288985), and [Bray et al](https://www.nature.com/articles/nbt.3519).

Remember that reads are broken into k-mers (here, k-mers of length 31). These k-mers are stored in a large and effective associative table, in which maximum mapable prefix are computed for each k-mer. Then, each probable position if compared in order to find a longest common prefix.

In short: reads are never mapped, but k-mer built from the reads are searched in a large associative table.

## What is the Expected-Maximization Algorithm?

[Expected maximization](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) algorithm is a method used to estimate the transcript abundances based on previously observed data within the sample. The main idea of this approach stands in two points:

    For each read, given a range of positions where this read comes from, guess the probability this read comes from a given position
    Knowing previous assignations, adjust the guessing-pattern

It has been shown by [Capellini et al.](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-1809.1955.tb01360.x), then by [Dempster et al.](https://www.jstor.org/stable/2984875), that EM algorithm converges when steps one and two are repeatedly applied one after each other.

## What is Variational Bayesian Expected-Maximisation Algorithm?

From [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#useem), we can find the following explanation about differences between EM algorithm, and the [VBEM](https://en.wikipedia.org/wiki/Variational_Bayesian_methods):

The details of the VBEM algorithm can be found in the article of [Patro et al.](https://www.nature.com/articles/nmeth.4197) While both the standard EM and the VBEM produce accurate abundance estimates, there are some trade-offs between the approaches. Specifically, the sparsity of the VBEM algorithm depends on the prior that is chosen. When the prior is small, the VBEM tends to produce a sparser solution than the EM algorithm, while when the prior is relatively larger, it tends to estimate more non-zero abundances than the EM algorithm. It is an active research effort to analyze and understand all the tradeoffs between these different optimization approaches. Also, the VBEM tends to converge after fewer iterations, so it may result in a shorter runtime; especially if you are computing many bootstrap samples.

## What are TPM? Why not use Raw couts, FPKM or RPM ?

TPM, FPKM, RPKM, and RPM are metrics which attempt to normalize the raw read count into a score. A very good [youtube video](https://www.youtube.com/watch?v=TTUrtCY2k-w) explains it and the following text repeats the video’s message:

In order to compare two values, the have to refer to the same scale. It is the very same process when you try to compare to academic averages: (12/20 and 3/5 are equal grades).

Here, in RNA-Seq, we have similar issues. The main issues are:

1. There is not the very same amount of biological material in all sequenced samples
1. Longer genes take longer time to be replicated, they will have fewer reads

RPM (Reads per Milion), only solved the first issue described above. To compute your RPM, you have to:

1. Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
1. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)

This metric has quickly been deprecated, yet many teams around the world still uses it, ignoring elementary PCR bias.

The idea was then to correct the read count by the gene length. There, the RPKM (Read Per Kilobase per Million) was born. Here’s how you do it for RPKM:

1. Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
1. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
1. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

Applied to Paired-end sequencing, this normalisation relies not on read-counts anymore, but on a pair-of-read counts. A pair of reads is called a fragment. There, the FPKM (Fragment Per Kilobase per Million) was born.

With RPKM or FPKM, the sum of normalized reads in each sample can be different. It means FPKM and RPKM does not solve the problem. They just reduce it.

Finally, TPM (Transcript Per Million), just switches two steps of the FPKM method, but solves the two issues presented at the beginning of this section. Here’s how you calculate TPM:

1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
1. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
1. Divide the RPK values by the “per million” scaling factor. This gives you TPM.

[Kanitz et al.](file:///home/tdayris/Documents/Projects/B17_FROL/Treatment_modif/report.html#NormalArticle1), and [Dillies et al.](file:///home/tdayris/Documents/Projects/B17_FROL/Treatment_modif/report.html#NormalArticle2) compares various methods of normalization including RPKM and TPM. If you persist in willing to use FPKM over TPM, you are doing things wrong when you compare two samples. Yet, you can use FPKM to compare two genes in one sample.

## What is the argument `--gencode` in Salmon?

This is a purely aesthetic option. Usually a genic object is identified with either: a gene name from [HUGO](https://www.genenames.org/) (i.e. PTEN, BRCA1, …), or from [ENSEMBL](https://www.ensembl.org/index.html) (i.e. ENSG00000171862, ENSG00000012048, …), etc. In its transcript genome sequences, [Gencode](https://www.gencodegenes.org/) uses multiple identifiers and separates them with a pipe (|).

As this is very convenient to search for names, transcripts, and protein correspondence, this also breaks compatibility with downstream tools. Salmon automatically breaks Gencodes names into a unique Ensembl identifier.
