#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This script contains functions that are to be called by any other scripts in
this pipeline.
"""

import argparse  # Argument parsing
import logging  # Logging behaviour
import pandas  # Handle large datasets
import pytest
import yaml  # Handle Yaml IO

import os.path as op  # Path and file system manipulation
import pandas  # Deal with TSV files (design)

from itertools import chain  # Chain iterators
from pathlib import Path  # Easily handle paths
from typing import Any, Dict, List, Optional, Union  # Type hints


# Building custom class for help formatter
class CustomFormatter(
        argparse.RawDescriptionHelpFormatter,
        argparse.ArgumentDefaultsHelpFormatter
    ):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """


def write_yaml(output_yaml: Path, data: Dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file
    """
    with output_yaml.open("w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)


def read_aggregation_table(count: str) -> pandas.DataFrame:
    """
    Load an aggregation table and return a DataFrame
    """
    # Build io paths objects
    count_table = Path(count)

    # Load dataset
    data = pandas.read_csv(
        count_table, sep="\t",
        index_col=0,
        header=0
    )

    logging.debug("Head of the count data:")
    logging.debug(data.head())

    return data


def fq_link(design: pandas.DataFrame) -> Dict[str, str]:
    """
    This function takes the "samples" described in config and returns
    a dictionnary with:
    sample file name : sample path
    """
    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired-ended case
        fq_list = chain(design["Upstream_file"], design["Downstream_file"])
    except KeyError:
        # Single ended case
        fq_list = design["Upstream_file"]

    return {op.basename(fq): op.realpath(fq) for fq in fq_list}


@pytest.mark.parametrize(
    "design, expected",
    [
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Upstream_file": "/path/to/S1_R1.fq",
                        "Downstream_file": "/other/to/S1_R2.fq",
                    }
                }
            ).T,
            {"S1_R1.fq": "/path/to/S1_R1.fq", "S1_R2.fq": "/other/to/S1_R2.fq"},
        ),
        (
            pandas.DataFrame({"S1": {"Upstream_file": "/path/to/S1_R1.fq"}}).T,
            {"S1_R1.fq": "/path/to/S1_R1.fq"},
        ),
        (
            pandas.DataFrame({"S1": {"Upstream_file": "path/to/S1_R1.fq"}}).T,
            {"S1_R1.fq": op.abspath("path/to/S1_R1.fq")},
        ),
    ],
)
def test_fq_link(
        design: pandas.DataFrame,
        expected: Optional[Dict[str, str]]
    ) -> None:
    """
    Test the function fq_link with multiple data
    """
    assert fq_link(design) == expected


def fq_root(design: pandas.DataFrame) -> Dict[str, str]:
    """
    This function takes the fastq file list and returns the root
    name corresponding to a fastq file
    sample name: sample link path
    """
    # For now, bz2 compression is not taken into account.
    possible_ext = ("fq", "fastq", "fq.gz", "fastq.gz")

    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired-ended case
        fq_list = chain(design["Upstream_file"], design["Downstream_file"])
    except KeyError:
        # Single ended case
        fq_list = design["Upstream_file"]

    # Build final result
    result = {}
    for fq_file in fq_list:
        # I always love writing these crazy for-break-else!
        for ext in possible_ext:
            if fq_file.endswith(ext):
                # Extension removal
                base = op.basename(fq_file)[: -(len(ext) + 1)]
                result[base] = f"raw_data/{op.basename(fq_file)}"
                break
        else:
            raise ValueError(f"Could not remove ext: {fq_file}")

    return result


@pytest.mark.parametrize(
    "design, expected",
    [
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Upstream_file": "/path/to/S1_R1.fq.gz",
                        "Downstream_file": "/other/to/S1_R2.fq.gz",
                    }
                }
            ).T,
            {"S1_R1": "raw_data/S1_R1.fq.gz", "S1_R2": "raw_data/S1_R2.fq.gz"},
        ),
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Upstream_file": "/path/to/S1_R1.fastq.gz",
                        "Downstream_file": "/other/to/S1_R2.fastq.gz",
                    }
                }
            ).T,
            {
                "S1_R1": "raw_data/S1_R1.fastq.gz",
                "S1_R2": "raw_data/S1_R2.fastq.gz",
            },
        ),
        (
            pandas.DataFrame(
                {"S1": {"Upstream_file": "/path/to/S1_R1.fastq.gz"}}
            ).T,
            {"S1_R1": "raw_data/S1_R1.fastq.gz"},
        ),
        (
            pandas.DataFrame(
                {"S1": {"Upstream_file": "/path/to/S1_R1.fastq"}}
            ).T,
            {"S1_R1": "raw_data/S1_R1.fastq"},
        ),
    ],
)
def test_fq_root(design: pandas.DataFrame, expected: Dict[str, str]) -> None:
    """
    Test the function fq_root with multiple data
    """
    assert fq_root(design) == expected


def ref_link(config: Dict[str, Any]) -> Dict[str, str]:
    """
    This function takes the "ref" described in config and returns
    a dictionnary with:
    ref file name : ref path
    """
    fasta = config["ref"]["fasta"]
    gtf = config["ref"]["gtf"]
    return {
        op.basename(fasta): op.realpath(fasta),
        op.basename(gtf): op.realpath(gtf),
    }


@pytest.mark.parametrize(
    "config, expected",
    [
        (
            {"ref": {"fasta": "/path/to/fasta.fa", "gtf": "/path/to/gtf.gtf"}},
            {"fasta.fa": "/path/to/fasta.fa", "gtf.gtf": "/path/to/gtf.gtf"},
        ),
        (
            {"ref": {"fasta": "path/to/fasta.fa", "gtf": "path/to/gtf.gtf"}},
            {
                "fasta.fa": op.abspath("path/to/fasta.fa"),
                "gtf.gtf": op.abspath("path/to/gtf.gtf"),
            },
        ),
    ],
)
def test_ref_link(config: Dict[str, Any], expected: Dict[str, str]) -> None:
    """
    Test the function ref_link with multiple arguments
    """
    assert ref_link(config) == expected


def fq_pairs(design: pandas.DataFrame) -> Dict[str, str]:
    """
    This function returns a sample ID and
    the corresponding fastq files.
    """
    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired end case
        iterator = zip(
            design["Sample_id"],
            design["Upstream_file"],
            design["Downstream_file"],
        )
        return {
            name: {
                "r1": f"raw_data/{op.basename(fq1)}",
                "r2": f"raw_data/{op.basename(fq2)}",
            }
            for name, fq1, fq2 in iterator
        }
    except KeyError:
        # Single end case
        iterator = zip(design["Sample_id"], design["Upstream_file"])
        return {
            name: {"r": f"raw_data/{op.basename(fq1)}"}
            for name, fq1 in iterator
        }


@pytest.mark.parametrize(
    "design, expected",
    [
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Sample_id": "S1",
                        "Upstream_file": "/path/to/S1_R1.fq.gz",
                        "Downstream_file": "/other/to/S1_R2.fq.gz",
                    }
                }
            ).T,
            {
                "S1": {
                    "r1": "raw_data/S1_R1.fq.gz",
                    "r2": "raw_data/S1_R2.fq.gz",
                }
            },
        ),
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Sample_id": "S1",
                        "Upstream_file": "/path/to/S1_R1.fq.gz",
                    }
                }
            ).T,
            {"S1": {"r": "raw_data/S1_R1.fq.gz"}},
        ),
        (
            pandas.DataFrame(
                {
                    "S1": {
                        "Sample_id": "S1",
                        "Upstream_file": "/path/to/S1_R1.fq.gz",
                        "Downstream_file": "/other/to/S1_R2.fq.gz",
                    },
                    "S2": {
                        "Sample_id": "S2",
                        "Upstream_file": "/path/to/S2_R1.fq.gz",
                        "Downstream_file": "/other/to/S2_R2.fq.gz",
                    },
                }
            ).T,
            {
                "S1": {
                    "r1": "raw_data/S1_R1.fq.gz",
                    "r2": "raw_data/S1_R2.fq.gz",
                },
                "S2": {
                    "r1": "raw_data/S2_R1.fq.gz",
                    "r2": "raw_data/S2_R2.fq.gz",
                },
            },
        ),
    ],
)
def test_fq_pairs(design: pandas.DataFrame, expected: Dict[str, Any]) -> None:
    """
    Test he function fq_pairs with multiple arguments
    """
    assert fq_pairs(design) == expected


def refs_pack(config: Dict[str, Any]) -> Dict[str, str]:
    """
    Return a dictionnary with references
    """
    return {
        "fasta": f"genomes/{op.basename(config['ref']['fasta'])}",
        "gtf": f"genomes/{op.basename(config['ref']['gtf'])}",
    }


@pytest.mark.parametrize(
    "config, expected",
    [
        (
            {"ref": {"fasta": "/path/to/fasta.fa", "gtf": "/path/to/gtf.gtf"}},
            {"fasta": "genomes/fasta.fa", "gtf": "genomes/gtf.gtf"},
        )
    ],
)
def test_ref_pack(config: Dict[str, Any], expected: Dict[str, str]) -> None:
    """
    Test the refs_pack function with multiple arguments
    """
    assert refs_pack(config) == expected


def salmon_quant_extra(config: Dict[str, Any]) -> str:
    """
    Return the corrected list of parameters for kallist quant
    """
    base = config["params"].get("salmon_quant_extra", "")
    return f"{base} --geneMap genomes/{op.basename(config['ref']['gtf'])}"


@pytest.mark.parametrize(
    "config, expected",
    [
        (
            {"params": {}, "ref": {"gtf": "genomes/gtf.gtf"}},
            " --geneMap genomes/gtf.gtf",
        ),
    ],
)
def test_salmon_quant_extra(config: Dict[str, Any], expected: str) -> None:
    """
    Test the salmon_quant_extra function with multiple arguments
    """
    assert salmon_quant_extra(config) == expected


def salmon_quant_output(config: Dict[str, Any]) -> str:
    """
    Return optionnal quant.genes if required by user
    """
    base = {"quant": "pseudo_mapping/{sample}/quant.sf"}

    try:
        # Case there is a GTF passed in the config file
        if config["ref"]["gtf"] != "" and config["ref"]["gtf"] is not None:
            base["quant_genes"] = "pseudo_mapping/{sample}/quant.genes.sf"
    except KeyError:
        pass

    return base


@pytest.mark.parametrize(
    "config, expected",
    [
        (
            {"ref": {"gtf": "genomes/annot.gtf"}},
            {
                "quant": "pseudo_mapping/{sample}/quant.sf",
                "quant_genes": "pseudo_mapping/{sample}/quant.genes.sf",
            },
        ),
        ({"ref": {}}, {"quant": "pseudo_mapping/{sample}/quant.sf"}),
        ({"ref": {"gtf": None}}, {"quant": "pseudo_mapping/{sample}/quant.sf"}),
    ],
)
def test_salmon_quant_output(config: Dict[str, Any], expected: str) -> None:
    """
    Test the salmon_quant_output function
    """
    assert salmon_quant_output(config) == expected


def sample_id(design: pandas.DataFrame) -> List[str]:
    """
    Return the list of samples identifiers
    """
    return design["Sample_id"].tolist()
