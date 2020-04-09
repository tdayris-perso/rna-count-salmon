"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""

import os               # OS related operations
import os.path as op    # Path and file system manipulation
import sys              # System related operations


from typing import Any, Dict, List     # Give IO information
import pandas as pd                    # Deal with TSV files (design)
from snakemake.utils import validate   # Check Yaml/TSV formats

try:
    from common_rna_count_salmon import *
except ImportError:
    print(locals())
    raise

# Snakemake-Wrappers version
swv = "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/0.51.0"
# github prefix
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/Unofficial"

# Loading configuration
if config == dict():
    configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pd.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")

report: "../report/general.rst"

def fq_pairs_w(wildcards) -> Dict[str, str]:
    """
    Dynamic wildcards call for snakemake.
    """
    try:
        return {**fq_pairs_dict[wildcards.sample],
                "gtf": refs_pack_dict["gtf"],
                "index": "salmon_index/genome_index"}
    except KeyError:
        return {**fq_pairs_dict[wildcards.sample],
                "index": "salmon_index/genome_index"}


def get_rcs_targets(get_fastqc: bool = False,
                    get_aggreg: bool = False,
                    get_multiqc: bool = False,
                    get_renamed: bool = False,
                    get_quant: bool = False) -> Dict[str, Any]:
    """
    This function returns the targets of Snakemake
    following the requests from the user.
    """
    targets = {}
    if config["workflow"]["fastqc"] is True and get_fastqc is True:
        targets["fastqc"] = expand(
            "qc/fastqc/{samples}_fastqc.{ext}",
            samples=fq_root_dict.keys(),
            ext=["html", "zip"]
        )

    if config["workflow"]["multiqc"] is True and get_multiqc is True:
        targets["multiqc"] = "qc/multiqc_report.html"

    if config["workflow"]["aggregate"] is True and get_aggreg is True:
        targets["aggregation"] = expand(
            "pseudo_mapping/aggregation/TPM{ext}.tsv",
            ext=["", ".genes"]
        )

    if get_renamed is True:
        targets["quant_renamed"] = expand(
            "pseudo_mapping/{sample}/quant.{sample}.tsv",
            sample=sample_id_list
        )

    if get_quant is True:
        targets["quant"] = expand(
            "pseudo_mapping/{sample}/quant.sf",
            sample=sample_id_list
        )

    return targets


# We will use these functions multiple times. On large input datasets,
# pre-computing all of these makes Snakemake faster.
fq_link_dict = fq_link(design)
fq_root_dict = fq_root(design)
ref_link_dict = ref_link(config)
fq_pairs_dict = fq_pairs(design)
refs_pack_dict = refs_pack(config)
sample_id_list = sample_id(design)
targets_dict = get_rcs_targets()

# Wildcards constraints
sample_constraint = "|".join(sample_id_list)
gene_ext_constraint = "|".join(["sf", "genes.sf"])
