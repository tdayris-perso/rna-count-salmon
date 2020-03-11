#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script aggregates multiple custom_data from MultiQC compatible
yaml files, and returns a single configuration file for MultiQC.
"""

import logging  # Traces and loggings
import numpy    # Handling maths on vectors
import pandas   # Handling large datasets
import yaml     # Handling yaml IO
import sys      # Handle system interactions

from pathlib import Path              # Easily deal with paths
from snakemake.utils import makedirs  # Create directories recursively
from typing import Any, Dict          # Type hints

from common_script_rna_count_salmon import write_yaml         # Common operations in this pipeline


def load_qc_config(qc: Path) -> Dict[str, Any]:
    """
    Read a MultiQC config file, and returns a dictionnary
    that is to be merged laterly
    """
    with qc.open("r") as config:
        return yaml.load(config, yaml.SafeLoader)


def multiqc_config(merged: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function builds a MultiQC-compatible dictionnary
    this header shall be added
    """
    return {
        "title": "MultiQC report - rna-count-salmon pipeline",
        "subtitle": "Bioinformatics Platform, Gustave Roussy Institute",
        "intro_text": ("This pipeline and all details are available at the"
                       " <a href='https://github.com/tdayris-perso/rna"
                       "-count-salmon'>rna-count-salmon Pipeline "
                       "page</a>"),
        "report_comment": "Powered by Snakemake",
        "custom_data": merged
    }


def test_multiqc_config() -> None:
    """
    Test the dict to multiqc convertion
    """
    test = {"custom_data": "ok"}
    expected = {
        "title": "MultiQC report - rna-count-salmon pipeline",
        "subtitle": "Bioinformatics Platform, Gustave Roussy Institute",
        "intro_text": ("MultiQC report summarise analysis results, "
                       "alongside with quality controls"),
        "report_comment": "Powered by ",
        "custom_data": "ok"
    }
    assert multiqc_config(test) == expected


def merge_configs(*configs: Dict[str, Any]) -> Dict[str, Any]:
    """
    This function merges multiple multiqc configurations and returns
    only one dictionnary
    """
    result = {}
    for section in configs:
        result.update(**section["custom_data"])
    logging.debug(result)

    return result


# def test_merge_configs

# The main process of the programm, not launched on import
if __name__ == '__main__':
    # Define logging behaviour
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )

    try:
        # Build IO paths
        configs = [Path(i) for i in snakemake.input["yaml"]]
        output = Path(snakemake.output["yaml"])

        # Handling yamls and merging them
        configs = [load_qc_config(config) for config in configs]
        merged = merge_configs(*configs)
        qc = multiqc_config(merged)
        logging.debug("Merged MultiQC data:")
        logging.debug(str(qc))

        # Output
        makedirs(str(output.parent))
        write_yaml(output, qc)
    except Exception as e:
        logging.exception("%s", e)
        sys.exit(1)

    sys.exit(0)
