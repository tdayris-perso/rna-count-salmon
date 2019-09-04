#!/usr/bin/python3.7
# conding: utf-8

"""
This script iterates over given salmon quantification files (quant.sf, or
quant.genes.sf). It extracts the columns "NumReads" in one hand, and "TPM"
in the other hand and creates two TSV-formatted text file with:

Column 1            - Target ID (transcript identifier or gene identifier)
Column 2 and others - Quantification for the given sample

First line of this file contains column header.

You can test this script with:
pytest -v aggregata_samples.py

This script is not intended to be used outside Snakemake.
"""

import logging            # Traces and loggings
import logging.handlers   # Logging behaviour
import os                 # Dealing with OS related function
import os.path as op      # Dealing with paths
import pandas as pd       # Parsing large quantification tables
import sys                # System related functions

from pathlib import Path              # Easily deal with paths
from snakemake.utils import makedirs  # Build output directories
from typing import List               # Type hints


logger = logging.getLogger(
    os.path.splitext(os.path.basename(sys.argv[0]))[0]
)


# Handling logging options
# No tests for this function
def setup_logging() -> None:
    """
    Configure logging behaviour
    """
    root = logging.getLogger("")
    root.setLevel(logging.WARNING)
    logger.setLevel(logging.DEBUG)
    if not args.quiet:
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(
            "%(levelname)s [%(name)s]: %(message)s"
        ))
        root.addHandler(ch)


# Parse quantification table
def read_salmon(path: str) -> pd.DataFrame:
    """
    Return a Pandas DataFrame from a salmon quantification file path

    Parameters:
        path    Path        Path to a salmon quantification result

    Return:
                DataFrame   A DataFrame containing the whole file

    Example:
    >>> read_salmon(Path("../tests/salmon_example/quant.sf"))
    """
    return pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float
        }
    )


def test_read_salmon() -> None:
    """
    This function tests Salmon reader function

    Example:
    pytest -v aggregata_samples.py -k test_read_salmon
    """
    test = read_salmon(Path("../tests/salmon_example/quant.2.sf")).sort_index()
    expected = pd.DataFrame({
        'Length': {
            'ENST00000387460.2': 66.0,
            'ENST00000387459.1': 69.0
        },
        'EffectiveLength': {
            'ENST00000387460.2': 23.0,
            'ENST00000387459.1': 24.0
        },
        'TPM': {
            'ENST00000387460.2': 30.3549,
            'ENST00000387459.1': 22.1085
        },
        'NumReads': {
            'ENST00000387460.2': 25.0,
            'ENST00000387459.1': 19.0
        }
    }).sort_index()

    print(test)
    print(expected)

    assert all(test == expected)


def extract_field(frame: pd.DataFrame,
                  sample_id: str,
                  column: str = "TPM") -> pd.DataFrame:
    """
    This function returns a single-columned dataframe from a full sample
    quantification dataframe.

    Parameters:
        frame       pd.DataFrame    The complete quantification dataframe
        sample_id   str             The sample id to be put as column name
        column      str             Then name of the column that is to be kept

    Return:
                    pd.DataFrame   The reduced dataframe, with only the
                                    selected column and sample renamed

    Example:
    >>> extract_field(example_frame, "id")
                            id
    Name
    ENST00000387460.2  30.3549
    ENST00000387459.1  22.1085
    """
    result = frame[[column]]
    result.columns = [sample_id]
    return result


def test_extract_field() -> None:
    """
    This function tests the extract_field function

    Example:
    pytest -v aggregata_samples.py -k test_extract_field
    """
    example_frame = pd.DataFrame({
        'EffectiveLength': {
            'ENST00000387459.1': 24.0,
            'ENST00000387460.2': 23.0
        },
        'Length': {
            'ENST00000387459.1': 69.0,
            'ENST00000387460.2': 66.0
        },
        'NumReads': {
            'ENST00000387459.1': 19.0,
            'ENST00000387460.2': 25.0
        },
        'TPM': {
            'ENST00000387459.1': 22.1085,
            'ENST00000387460.2': 30.3549
        }
    }).sort_index()

    expected = pd.DataFrame(
        {'id': {'ENST00000387459.1': 22.1085, 'ENST00000387460.2': 30.3549}}
    ).sort_index()

    assert all(expected == extract_field(example_frame, "id"))


def merge_reduced_frames(*paths: List[str],
                         prefix: str = "",
                         column: str = "TPM") -> pd.DataFrame:
    """
    Return a dataframe containing the required field

    Parameters:
        paths       List[str]           A List of Paths to each of the salmon
                                        quantification files
        prefix      str                 The prefix to remove from sample path
                                        in order to get sample id
        column      str                 The column that is to be kept

    Return:
                    pd.DataFrame        The merged DataFrame

    Example:
    >>> merge_reduced_frames("quant.sf", "quant.2.sf")
                       quant.sf  quant.2.sf
    Name
    ENST00000387460.2   16.4862     30.3549
    ENST00000387459.1   24.4589     22.1085
    """
    merged_frame = None
    for path in paths:
        logger.debug("Working on {}".format(path))
        data = extract_field(
            read_salmon(path),
            str(path) if prefix == "" else str(path)[len(prefix):],
            column
        )

        try:
            merged_frame = pd.merge(
                merged_frame,
                data,
                left_index=True,
                right_index=True
            )
        except TypeError:
            merged_frame = data
        except ValueError:
            merged_frame = data

    merged_frame.fillna(0)
    logger.debug("Head of the {}-merged frame:".format(column))
    logger.debug(merged_frame.head())
    return merged_frame


def test_merge_reduced_frames() -> None:
    """
    This function tests the merge_reduced_frames frunction

    Example:
    pytest -v aggregated_samples.py -k test_merge_reduced_frames
    """
    expected = pd.DataFrame(
         {'quant.sf': {
            'ENST00000387459.1': 24.4589,
            'ENST00000387460.2': 16.4862
         },
          'quant.2.sf': {
            'ENST00000387459.1': 22.1085,
            'ENST00000387460.2': 30.3549
          }}
    ).sort_index()

    got = merge_reduced_frames(
        "../tests/salmon_example/quant.sf",
        "../tests/salmon_example/quant.2.sf",
        prefix="../tests/salmon_example/"
    ).sort_index()

    assert all(got == expected)


# The main process of the programm, not launched on import
if __name__ == '__main__':
    # Building output directory
    makedirs(op.dirname(snakemake.output["NumReads"]))

    # Iterating through columns
    for column in ["NumReads", "TPM"]:
        data = merge_reduced_frames(
            *snakemake.input["quants"],
            prefix="pseudo_mapping/",
            column=column
        )

        logger.debug(data.head())

        # Saving to TSV formatted text file
        data.to_csv(snakemake.output[column], sep="\t")
