#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script extracts statistics on null counts based on aggregated
salmon's quant files (either genes or transcripts). It produces:

- A yaml file that suits MultiQC requirements
- A png image (histograms)
"""

# Maths
import pandas   # Handle large datasets
import numpy    # Handle vectors and maths

# Text IO
import logging  # Log activity
import yaml     # Handle Yaml IO
import sys      # System interactions

# Plotting
from bokeh.io import export_png     # Handle graphics IO
from bokeh.plotting import figure   # Handle plotting object
from bokeh.models import ColumnDataSource, HoverTool   # Optional behaviours

# Language
from pandas.testing import assert_frame_equal  # Dedicated assertion
from pathlib import Path              # Easily handle paths
from snakemake.utils import makedirs  # Make directories recursively
from typing import Any, Dict          # Type hints

# Common functions within pipeline
from common import write_yaml, read_aggregation_table


def per_sample_null_count(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    This function computes per sample statistics on null count
    and returns these results as a DataFrame

        Parameters:
            data    DataFrame       The count table

        Return:
                    DataFrame       Null counts proportions
    """
    # Count null number
    nb_null = (data == 0).sum().to_dict()
    nb_not_null = (data > 0).sum().to_dict()

    # Get null proportion
    nb_obs = len(data)
    prop_null = {
        k: v / nb_obs
        for k, v in nb_null.items()
    }

    null_counts = pandas.DataFrame({
        "nb_null": nb_null,
        "nb_not_null": nb_not_null,
        "prop_null": prop_null
    })

    logging.debug("Null counts:")
    logging.debug(null_counts.head())

    return null_counts


def test_per_sample_null_count() -> None:
    """
    This function tests the per_sample_null_count function.
    """
    frame = pandas.DataFrame({'a': [1, 2, 0], 'b': [0, 1, 0]})
    expected = pandas.DataFrame({
        'nb_null': {'a': 1, 'b': 2},
        'nb_not_null': {'a': 2, 'b': 1},
        'prop_null': {'a': 0.3333333333333333, 'b': 0.6666666666666666}
    })
    assert_frame_equal(per_sample_null_count(frame), expected)


def to_multiqc_yaml(data: pandas.DataFrame) -> Dict[str, Any]:
    """
    This function takes a null count statistics and returns a
    yaml that suits multiqc.
    """
    return {
        'custom_data': {
            'null_count_stats': {
                'id': 'null_count_section',
                'section_name': 'Null count statistics',
                'description': 'This graph shows null count per sample',
                'plot_type': 'bargraph',
                'pconfig': {
                    'id': 'null_count_histogram',
                    'title': 'Null count histogram',
                    'ylab': 'Number of null counts',
                    'xlab': 'Samples',
                },
                'data': data[["nb_null", "nb_not_null"]].to_dict()
            }
        }
    }


def test_to_multiqc_yaml() -> None:
    """
    This function tests the pandas to multiqc converter
    """
    frame = pandas.DataFrame({
        'nb_null': {'a': 1, 'b': 2},
        'nb_not_null': {'a': 2, 'b': 1},
        'prop_null': {'a': 0.3333333333333333, 'b': 0.6666666666666666}
    })
    expected = {
        'custom_data': {
            'null_count_stats': {
                'id': 'null_count_section',
                'section_name': 'Null count statistics',
                'description': 'This graph shows null count per sample',
                'plot_type': 'bargraph',
                'pconfig': {
                    'id': 'null_count_histogram',
                    'title': 'Null count histogram',
                    'ylab': 'Number of null counts',
                    'xlab': 'Samples'
                },
                'data': {
                    'nb_null': {'a': 1, 'b': 2},
                    'nb_not_null': {'a': 2, 'b': 1}
                }
            }
        }
    }
    assert to_multiqc_yaml(frame) == expected


def bokeh_figure(data: pandas.DataFrame) -> figure:
    """
    This function produces a bokeh.models.figure object
    which represents a histogram of null count data
    """
    # Reducing data to columns of interest
    cols = ["nb_null", "nb_not_null"]
    data = data[cols].reset_index()
    data.columns = ["sample"] + cols

    # Building data structures for bokeh
    samples = data.index.tolist()
    source = ColumnDataSource(data=data)

    hover = HoverTool(
        tooltips=[
            ("Sample", "@sample"),
            ("Null counts", "@nb_null"),
            ("Non-null counts", "@nb_not_null")
        ]
    )

    # Building figure
    p = figure(
        y_range=data.index.tolist(),
        plot_height=720,
        plot_width=1024,
        title="Per sample null count",
    )
    p.add_tools(hover)

    p.hbar_stack(
        cols,
        y="index",
        height=0.9,
        color=['#3182bd', '#9ecae1'],
        source=source
    )
    return p


def write_png(fig: figure, output_png: Path) -> None:
    """
    This function saves the bokeh plot as a png file
    """
    export_png(fig, filename=str(output_png))


if __name__ == '__main__':
    # Build logging object and behaviour
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )

    try:
        # Extract statistics
        data = read_aggregation_table(snakemake.input["count"])
        null_data = per_sample_null_count(data)
    except Exception as e:
        logging.exception("%s", e)
        sys.exit(1)

    # Output formatted data on demand
    try:
        output_yaml = Path(snakemake.output["yaml"])
        logging.debug("Yaml output:")
        logging.debug(output_yaml)
        write_yaml(output_yaml, to_multiqc_yaml(null_data))
    except KeyError as ke:
        logging.exception("%s", ke)

    try:
        fig = bokeh_figure(null_data)
    except Exception as e:
        logging.exception("%s", e)
        sys.exit(1)

    try:
        output_png = Path(snakemake.output["png"])
        makedirs(str(output_png.parent))
        write_png(fig, output_png)
    except KeyError as ke:
        logging.exception("%s", ke)
    sys.exit(0)
