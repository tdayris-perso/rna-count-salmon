import snakemake.utils  # Load snakemake API
import sys              # System related operations

# Python 3.8 is required
if sys.version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later.")

# Snakemake 5.14.0 at least is required
snakemake.utils.min_version("5.14.0")

include: "rules/common.smk"
include: "rules/copy.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/salmon.smk"
include: "rules/aggregation.smk"

workdir: config["workdir"]
container: config["singularity_docker_image"]
localrules: copy_fastq, copy_extra

rule all:
    input:
        **get_rcs_targets(get_multiqc=True, get_aggreg=True,
                          get_renamed=True, get_fastqc=True)
    message:
        "Finishing the Salmon RNA-Seq quantification pipeline"
