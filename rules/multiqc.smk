"""
This rule runs MultiQC in order to collect metrics on most of our tools and
raw files: Fastq + Salmon. We need to include the fasta reference for
the report option only.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html
"""
rule multiqc:
    input:
        **get_rcs_targets(get_fastqc=True, get_quant=True, get_qc_config=True)
    output:
        report(
            "qc/multiqc_report.html",
            caption="../report/multiqc.rst",
            category="Quality Controls"
        )
    params: "-c qc/multiqc_configs/complete_multiqc_config.yaml"
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 2048, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 120)
        )
    log:
        "logs/multiqc.log"
    message:
        "Gathering quality reports with MultiQC"
    wrapper:
        f"{swv}/bio/multiqc"


"""
This rule aggregates separate optional graphs in yaml format, and
builds the complete MultiQC configuration file
"""
rule aggregate_multiqc_config:
    input:
        yaml = [
            "qc/multiqc_configs/multiqc_null_count.yaml",
            "qc/multiqc_configs/multiqc_box_count.yaml"
        ]
    output:
        yaml = "qc/multiqc_configs/complete_multiqc_config.yaml"
    message:
        "Aggregating multiqc configuration"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(
                attempt * len(sample_id_list) * 250, 10240
            )
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 60)
        )
    log:
        "logs/multiqc/aggregate_config.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/multiqc_config.py"
