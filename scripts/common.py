#!/usr/bin/python3.7
# -*- coding: utf-8 -*-


"""
This script contains functions that are to be called by any other scripts in
this pipeline.
"""

import argparse        # Argument parsing
import logging         # Logging behaviour
import pandas          # Handle large datasets
import yaml            # Handle Yaml IO

from pathlib import Path             # Easily handle paths
from typing import Any, Dict, Union  # Type hints


# Building custom class for help formatter
class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """
    pass


def write_yaml(output_yaml: Path, data: Dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file
    """
    with output_yaml.open("w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)


def read_aggregation_table(count: str) -> pandas.DataFrame:
    # Build io paths objects
    count_table = Path(count)

    # Load dataset
    data = pandas.read_csv(
        count_table,
        sep="\t",
        index_col=0,
        header=0
    )

    logging.debug("Head of the count data:")
    logging.debug(data.head())

    return data
