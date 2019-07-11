"""
This rule computes the salmon index from a fasta-formatted
transcriptome sequence (no genome)

more information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/salmon/index.html
"""
rule salmon_index:
    input:
        **refs_pack_dict
    output:
        index = directory("salmon_index/genome_index")
    message:
        "Indexing {input.fasta} with Salmon"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 2048 + 10240, 35580)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 15 + 15, 180)
        )
    version: "1.0"
    threads:
        min(config["threads"], 12)
    params:
        extra = config["params"]["salmon_index_extra"]
    log:
        "logs/salmon/index.log"
    wrapper:
        f"{swv}/bio/salmon/index"

"""
This rule performs the quantification step with pseudo-mapping
fromfastq-formatted rnaseq reads and salmon genome index.

More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/salmon/quant.html
"""
rule salmon_quant:
    input:
        unpack(fq_pairs_w)
    output:
        quant = report(
            "pseudo_mapping/{sample}/quant.sf",
            caption="../report/salmon.transcripts.rst",
            category="Sample Count"
        )
    message:
        "Quantifying {wildcards.sample} with Salmon"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 5120 + 2048, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 15 + 15, 180)
        )
    version: "1.0"
    threads:
        min(config["threads"], 12)
    params:
        libType = config["params"]["libType"],
        extra = salmon_quant_extra()
    log:
        "logs/salmon/quant_{sample}.log"
    wrapper:
        f"{swv}/bio/salmon/quant"
