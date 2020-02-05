#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script takes a TSV-formatted count table and another TSV-formatted
correspondancy between genes and transcripts identifiers.
"""

import logging  # Traces and loggings
import pandas   # Handling large datasets
import sys      # Handle system interactions

from pathlib import Path                   # Easily deal with paths
from snakemake.utils import makedirs       # Create directories recursively
from common import read_aggregation_table  # Common operations in this pipeline


if __name__ == '__main__':
    # Define logging behaviour
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )

    try:
        # Define IO
        counts = Path(snakemake.input["counts"])
        tr2gene = Path(snakemake.input["tr2gene"])
        output = Path(snakemake.output["annotated"])

        # Load datasets
        data = read_aggregation_table(counts)
        logging.debug("Head of counts table:")
        logging.debug(data.head())
        annot = pandas.read_csv(
            tr2gene,
            sep="\t",
            header=0,
            index_col=None,
            dtype=str
        )
        logging.debug("Head of annotation table:")
        logging.debug(annot.head())

        # Annotate
        right_key = (
            "gene_id"
            if snakemake.params["on_genes"] is True
            else "transcript_id"
        )

        logging.debug(f"Right key used for annotation merge: {right_key}")
        result = pandas.merge(
            data,
            annot,
            left_on="Name",
            right_on=right_key
        )
        logging.debug("Head of merged table:")
        logging.debug(result.head())

        # Save results
        makedirs(str(output.parent))
        result.to_csv(
            output,
            sep="\t"
        )
    except Exception as e:
        logging.exception("%s", e)
        sys.exit(1)

    sys.exit(0)
