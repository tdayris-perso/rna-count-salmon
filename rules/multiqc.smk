"""
This rule runs MultiQC in order to collect metrics on most of our tools and
raw files: Fastq + Salmon. We need to include the fasta reference for
the report option only.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html
"""
rule multiqc:
    input:
        **get_rcs_targets(get_fastqc=True, get_quant=True)
    output:
        report(
            "qc/multiqc_report.html",
            caption="../report/multiqc.rst",
            category="Quality Controls"
        )
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 2048
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 45
        )
    log:
        "logs/multiqc.log"
    message:
        "Gathering quality reports with MultiQC"
    wrapper:
        f"{swv}/bio/multiqc"
