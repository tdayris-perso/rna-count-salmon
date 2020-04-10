"""
This rule computes the salmon index from a fasta-formatted
transcriptome sequence (no genome)

more information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/salmon/index.html
"""
rule salmon_index:
    input:
        fasta = refs_pack_dict['fasta']
    output:
        index = directory("salmon_index/genome_index")
    message:
        "Indexing {input.fasta} with Salmon"
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 2048 + 10240
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 15 + 60
        )
    threads:
        min(config["threads"], 12)
    params:
        extra = config["params"].get("salmon_index_extra", "")
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
        quant = "pseudo_mapping/{sample}/quant.sf",
        quant_genes = "pseudo_mapping/{sample}/quant.genes.sf"
    message:
        "Quantifying {wildcards.sample} with Salmon"
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 5120 + 2048
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 15 + 75
        )
    wildcard_constraints:
        sample = sample_constraint
    threads:
        min(config["threads"], 12)
    params:
        libType = config["params"].get("libType", "A"),
        extra = salmon_quant_extra(config)
    log:
        "logs/salmon/quant_{sample}.log"
    wrapper:
        f"{swv}/bio/salmon/quant"


"""
This rule solves the issue: Duplicate explicit target name: "quant.sf"
within reporting
"""
rule salmon_quant_rename:
    input:
        "pseudo_mapping/{sample}/quant.sf"
    output:
        "pseudo_mapping/{sample}/quant.{sample}.tsv"
    message:
        "Symbolic link for quantification of {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: attempt * 256 + 128
        ),
        time_min = (
            lambda wildcards, attempt: attempt * 3 + 2
        )
    log:
        "logs/salmon/rename/{sample}.log"
    conda:
        "../envs/bash.yaml"
    params:
        pwd = os.getcwd()
    wildcard_constraints:
        sample = sample_constraint
    shell:
        "ln --verbose --symbolic {params.pwd}/{input} {params.pwd}/{output} > {log} 2>&1"
