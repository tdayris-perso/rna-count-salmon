"""
This rule performs quantification aggregations: from multiple quant.sf to
a single large quantification table for all sample
"""
rule aggregate_salmon_gene_counts:
    input:
        quant = expand(
            "pseudo_mapping/{sample}/quant{ext}.sf",
            sample=sample_id_list,
            allow_missing=True
        ),
        tx2gene = "pseudo_mapping/aggregation/transcript_to_gene.tsv"
    output:
        tsv = report(
            "pseudo_mapping/aggregation/TPM{ext}.tsv",
            category="Counts",
            caption="../report/aggregated_counts.rst"
        )
    message:
        "Aggregating all targets{wildcards.ext} abundancies"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * len(sample_id_list) * 250
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 10
        )
    params:
        drop_null = True,
        drop_na = True,
        header = True,
        genes = (
            lambda wildcards: True if wildcards.ext == ".genes" else False
        ),
        index_label = True
    wildcard_constraints:
        ext = "|.genes"
    log:
        "logs/aggregation/aggregate_{ext}.log"
    wrapper:
        f"{git}/bio/pandas/salmon"


"""
This rule creates the transcript to gene table from a GTF file. It's a bed-like
TSV-formatted text file.
"""
rule tr2gene:
    input:
        gtf = refs_pack_dict["gtf"]
    output:
        tsv = temp("pseudo_mapping/aggregation/transcript_to_gene.tsv")
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
    params:
        header = True,
        positions = True
    log:
        "logs/transcript_to_gene.log"
    wrapper:
        f"{git}/bio/tx_to_gene/gtf"
