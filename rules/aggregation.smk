"""
This rule performs quantification aggregations: from multiple quant.sf to
a single large quantification table for all sample
"""
rule aggregate_salmon_transcripts_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}/quant.sf",
            sample=sample_id_list
        )
    output:
        NumReads = report(
            "aggregated_salmon_counts/NumReads_transcripts.tsv",
            caption="../report/raw_counts.rst",
            category="Aggregated Counts"
        ),
        TPM = report(
            "aggregated_salmon_counts/TPM_transcripts.tsv",
            caption="../report/tpm.rst",
            category="Aggregated Counts"
        )
    message:
        "Aggregating all salmon abundancies"
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
        "logs/aggregation.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/aggregate_samples.py"

"""
This rule performs quantification aggregations: from multiple quant.sf to
a single large quantification table for all sample
"""
rule aggregate_salmon_gene_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}/quant.genes.sf",
            sample=sample_id_list
        )
    output:
        NumReads = report(
            "aggregated_salmon_counts/NumReads_genes.tsv",
            caption="../report/raw_counts.rst",
            category="Aggregated Counts"
        ),
        TPM = report(
            "aggregated_salmon_counts/TPM_genes.tsv",
            caption="../report/tpm.rst",
            category="Aggregated Counts"
        )
    message:
        "Aggregating all salmon abundancies"
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
        "logs/aggregation.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/aggregate_samples.py"
