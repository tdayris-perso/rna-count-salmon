"""
This rule collects metrics on raw fastq files.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html
"""
rule fastqc:
    input:
        lambda wildcards: fq_root_dict[wildcards.sample]
    output:
        html = report(
            "qc/fastqc/{sample}_fastqc.html",
            caption="../report/fastqc.rst",
            category="Quality Controls"
        ),
        zip = "qc/fastqc/{sample}_fastqc.zip"
    params:
        ""
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 8096)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 120)
        )
    log:
        "logs/fastqc/{sample}.log"
    message:
        "Controling quality of {wildcards.sample} fastq file with FastQC"
    wrapper:
        f"{git}/bio/fastqc"
