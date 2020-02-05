#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script extracts statistics on non-null counts based on aggregated
salmon's quant files (either genes or granscripts). It produces:

- A yaml file that suits MultiQC requirements
- A png image (box plots)
"""


# Maths
import pandas  # Handle large datasets
import numpy   # Handle vectors and maths

# Test IO
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
from typing import Any, Dict, List    # Type hints

# Common functions within pipeline
from common import write_yaml, read_aggregation_table


def per_sample_quantiles(data: pandas.DataFrame) -> pandas.DataFrame:
    """
    This function computes per sample statistics on non-null counts
    and returns these results as a DataFrame

    Parameters:
        data    DataFrame       The count table

    Return:
                DataFrame       Count quantiles
    """
    # Filter full-null lines
    not_null = data.loc[~(data == 0).all(axis=1)]

    # Compute quantiles
    quantiles = not_null.quantile(
        q=[0.25, 0.5, 0.75],
        axis=0,
        numeric_only=True,
        interpolation='linear'
    )

    mean = not_null.mean(
        axis=0
    )

    result = pandas.concat(
        [quantiles.T, mean],
        axis=1
    ).T

    result["title"] = [
        "Lower Quartile",
        "Median",
        "Higher Quartile",
        "Mean"
    ]

    result.set_index(
        "title",
        inplace=True
    )
    del result.index.name

    logging.debug(f"Statistics on non-null counts:")
    logging.debug(result.head())

    return result


def test_per_sample_quantiles() -> None:
    """
    This function tests the per_sample_quantiles function.
    """
    test = pandas.DataFrame(
        {'a': [1, 2, 0], 'b': [0, 1, 0]}
    )
    expected = pandas.DataFrame(
        {'a': {"Lower Quartile": 1.25,
               "Median": 1.50,
               "Higher Quartile": 1.75,
               "Mean": 1.50},
         'b': {"Lower Quartile": 0.25,
               "Median": 0.50,
               "Higher Quartile": 0.75,
               "Mean": 0.50}}
    )
    expected = expected.T[
        ['Lower Quartile', 'Median', 'Higher Quartile', 'Mean']
    ].T

    assert_frame_equal(per_sample_quantiles(test), expected)


def to_multiqc_yaml(data: pandas.DataFrame) -> Dict[str, Any]:
    """
    This function takes the pahdas.DataFrame and builds a multiqc
    compatible dictionnary.
    """
    return {
        'custom_data': {
            'my_data_type': {
                'id': 'quantiles_section',
                'section_name': 'Count quantiles',
                'description': 'These box plots shows counts across samples',
                'plot_type': 'table',
                'pconfig': {
                    'id': 'quantiles_section',
                    'namespace': 'Quantiles'
                },
                'headers': {
                    'Lower Quartile': {
                        "description": 'Value of the lower counts quartile',
                        "format": "{:,.2f}"
                    },
                    'Medians': {
                        "description": 'Value of the mean counts',
                        "format": "{:,.2f}"
                    },
                    'Higher Quartile': {
                        "description": 'Value of the higher counts quartile',
                        "format": "{:,.2f}"
                    },
                    'Mean': {
                        "description": "Value of the median count",
                        "format": "{:,.2f}"
                    }
                },
                'data': data.to_dict()
            }
        }
    }


def test_to_multiqc_yaml() -> None:
    """
    This function tests to_multiqc_yaml function
    """
    test = pandas.DataFrame(
        {'a': {"Lower Quartile": 1.25,
               "Median": 1.5,
               "Higher Quartile": 1.75,
               "Mean": 1.0},
         'b': {"Lower Quartile": 0.25,
               "Median": 0.5,
               "Higher Quartile": 0.75,
               "Mean": 0.0}}
    )
    expected = {
        'custom_data': {
            'box_count_stats': {
                'id': 'quantiles_section',
                'section_name': 'Count quantiles',
                'description': 'These box plots shows counts across samples',
                'plot_type': 'table',
                'pconfig': {
                    'id': 'quantiles_section',
                    'namespace': 'Quantiles'
                },
                'headers': {
                    'Lower Quartile': {
                        "description": 'Value of the lower counts quartile',
                        "format": "{:,.2f}"
                    },
                    'Medians': {
                        "description": 'Value of the mean counts',
                        "format": "{:,.2f}"
                    },
                    'Higher Quartile': {
                        "description": 'Value of the higher counts quartile',
                        "format": "{:,.2f}"
                    },
                    'Mean': {
                        "description": "Value of the median count",
                        "format": "{:,.2f}"
                    }
                },
                'data': {
                    'a': {
                        "Lower Quartile": 1.25,
                        "Median": 1.5,
                        "Higher Quartile": 1.75,
                        "Mean": 1.0},
                    'b': {
                        "Lower Quartile": 0.25,
                        "Median": 0.5,
                        "Higher Quartile": 0.75,
                        "Mean": 0.0
                    }
                }
            }
        }
    }
    assert to_multiqc_yaml(test) == expected


if __name__ == '__main__':
    # Build logging object and behaviour
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG
    )

    try:
        # Extract non-null counts and statistics
        data = read_aggregation_table(snakemake.input["count"])
        qtiles = per_sample_quantiles(data)
    except Exception as e:
        logging.exception("%s", e)
        sys.exit(1)

    # Output formatted data on demand
    try:
        output_yaml = Path(snakemake.output["yaml"])
        logging.debug("Yaml output:")
        logging.debug(output_yaml)
        write_yaml(output_yaml, to_multiqc_yaml(qtiles))
    except KeyError as ke:
        logging.exception("%s", ke)
    sys.exit(0)
