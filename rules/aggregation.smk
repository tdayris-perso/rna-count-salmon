"""
This rule performs quantification aggregations: from multiple quant.sf to
a single large quantification table for all sample
"""
rule aggregate_salmon_transcripts_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}/quant.{{ext}}",
            sample=sample_id_list
        )
    output:
        NumReads = temp("aggregated_salmon_counts/NumReads.{ext}.tsv"),
        TPM = temp("aggregated_salmon_counts/TPM.{ext}.tsv")
    message:
        "Aggregating all abundancies (on {wildcards.ext})"
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
        "logs/aggregation/aggregate_{ext}.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/aggregate_samples.py"


"""
This rule computes proportion of null counts over samples. This is a basic
quality control.
"""
rule null_counts:
    input:
        count = "aggregated_salmon_counts/NumReads.sf.tsv"
    output:
        yaml = "qc/multiqc_configs/multiqc_null_count.yaml",
        png = "qc/images/null_counts.png"
    message:
        "Controling quality over null counts proportion in all samples"
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
        "logs/null_counts.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/null_counts.py"


"""
This rule performs additional statistics on per-sample gene-count quartiles
"""
rule box_counts:
    input:
        count = "aggregated_salmon_counts/NumReads.sf.tsv"
    output:
        yaml = "qc/multiqc_configs/multiqc_box_count.yaml"
    message:
        "Controling quality over counts' median and quartiles in all samples"
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
        "logs/box_counts.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/box_counts.py"


"""
This rule creates the transcript to gene table from a GTF file. It's a bed-like
TSV-formatted text file.
"""
rule tr2gene:
    input:
        gtf = refs_pack_dict.get("gtf", "")
    output:
        tsv = temp("aggregated_salmon_counts/tr2gene.tsv")
    message:
        "Creating a transcript to gene table based on GTF"
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
        "logs/transcript_to_gene.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/transcript_to_gene.py"


"""
This rule adds genes information in the count tables based on the
transcript to gene table
"""
rule annotate_counts:
    input:
        counts = "aggregated_salmon_counts/{content}.{ext}.tsv",
        tr2gene = "aggregated_salmon_counts/tr2gene.tsv"
    output:
        annotated = report(
            "aggregated_salmon_counts/{content}.{ext}.annotated.tsv",
            caption="../report/aggregated_counts.rst",
            category="Aggregated Counts"
        )
    message:
        "Adding genes information in {wildcards.content} ({wildcards.ext})"
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
    wildcard_constraints:
        content = "TPM|NumReads",
        ext = "sf|genes.sf"
    params:
        on_genes = (
            lambda w: w.ext == "genes.sf"
        )
    log:
        "logs/annotate/{content}/{ext}.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/annotate_counts.py"
