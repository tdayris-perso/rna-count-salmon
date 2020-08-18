#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script aims to prepare the configuration file used
by the rna-count-salmon pipeline

It goes through the arguments passed in command line and
builds a yaml formatted text file used as a configuration
file for the snakemake pipeline.

You can test this script with:
pytest -v ./prepare_config.py

Usage example:
# Whole pipeline
python3.7 ./prepare_config.py /path/to/fasta_file.fa

# No quality controls, only quantification
python3.7 ./prepare_config.py /path/to/fasta_file.fa --no-fastqc --no-multiqc

# Whole pipeline, verbose mode activated
python3.7 ./prepare_config.py /path/to/fasta_file.fa -v
"""


import argparse  # Parse command line
import logging  # Traces and loggings
import os  # OS related activities
import pytest  # Unit testing
import shlex  # Lexical analysis
import sys  # System related methods
import yaml  # Parse Yaml files

from pathlib import Path  # Paths related methods
from snakemake.utils import makedirs  # Easily build directories
from typing import Dict, Any  # Typing hints

try:
    from scripts.common_script_rna_count_salmon import *
except ModuleNotFoundError:
    from common_script_rna_count_salmon import *


def parser() -> argparse.ArgumentParser:
    """
    Build the argument parser object
    """
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does not make any magic. Please check the prepared"
        " configuration file!",
    )

    # Parsing positional argument
    main_parser.add_argument(
        "fasta",
        help="Path to the fasta-formatted transcriptome sequence",
        type=str,
    )

    main_parser.add_argument(
        "gtf",
        help="Path to GTF-formatted genome annotation",
        type=str
    )

    # Parsing optional arguments
    main_parser.add_argument(
        "--design",
        help="Path to design file (default: %(default)s)",
        type=str,
        metavar="PATH",
        default="design.tsv",
    )

    main_parser.add_argument(
        "--workdir",
        help="Path to working directory (default: %(default)s)",
        type=str,
        metavar="PATH",
        default=".",
    )

    main_parser.add_argument(
        "--threads",
        help="Maximum number of threads used (default: %(default)s)",
        type=int,
        default=1,
    )

    main_parser.add_argument(
        "--singularity",
        help="Docker/Singularity image (default: %(default)s)",
        type=str,
        default="docker://continuumio/miniconda3:4.4.10",
    )

    main_parser.add_argument(
        "--cold-storage",
        help="Space separated list of absolute path to "
        "cold storage mount points (default: %(default)s)",
        nargs="+",
        type=str,
        default=[" "],
    )

    main_parser.add_argument(
        "--no-fastqc",
        help="Do not perform any fastqc",
        action="store_true"
    )

    main_parser.add_argument(
        "--no-multiqc",
        help="Do not perform final multiqc",
        action="store_true"
    )

    main_parser.add_argument(
        "--aggregate",
        help="Perform sample count aggregation",
        action="store_true",
    )

    main_parser.add_argument(
        "--salmon-index-extra",
        help="Extra parameters for salmon index step (default: %(default)s)",
        type=str,
        default="--keepDuplicates --gencode --perfectHash",
    )

    main_parser.add_argument(
        "--salmon-quant-extra",
        help="Extra parameters for salmon quantification step "
        "(default: %(default)s)",
        type=str,
        default="--numBootstraps 100 --validateMappings --gcBias --seqBias",
    )

    main_parser.add_argument(
        "--libType",
        help="The salmon library type (default: %(default)s)",
        type=str,
        default="A",
    )

    # Logging options
    log = main_parser.add_mutually_exclusive_group()
    log.add_argument(
        "-d",
        "--debug",
        help="Set logging in debug mode",
        default=False,
        action="store_true",
    )

    log.add_argument(
        "-q",
        "--quiet",
        help="Turn off logging behaviour",
        default=False,
        action="store_true",
    )

    return main_parser


# Argument parsing functions
def parse_args(args: Any) -> argparse.ArgumentParser:
    """
    This function parses command line arguments

    Parameters
        args     Any             All command line arguments

    Return
                ArgumentParser   A object designed to parse the command line

    Example:
    >>> parse_args(shlex.split("/path/to/fasta --no-fastqc"))
    Namespace(aggregate=False, cold_storage=[' '], debug=False,
    design='design.tsv', fasta='/path/to/fasta', gtf=None, libType='A',
    no_fastqc=False, no_multiqc=False, quiet=False, salmon_index_extra='
    --keepDuplicates --gencode --perfectHash', salmon_quant_extra='
    --numBootstraps 100 --validateMappings --gcBias --seqBias',
    singularity='docker://continuumio/miniconda3:4.4.10',
    threads=1, workdir='.')
    """
    return parser().parse_args(args)


def test_parse_args() -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    options = parse_args(shlex.split("/path/to/fasta /path/to/gtf"))
    expected = argparse.Namespace(
        aggregate=False,
        cold_storage=[" "],
        debug=False,
        design="design.tsv",
        fasta="/path/to/fasta",
        gtf="/path/to/gtf",
        libType="A",
        no_fastqc=False,
        no_multiqc=False,
        quiet=False,
        salmon_index_extra="--keepDuplicates --gencode --perfectHash",
        salmon_quant_extra=(
            "--numBootstraps 100 --validateMappings " "--gcBias --seqBias"
        ),
        singularity="docker://continuumio/miniconda3:4.4.10",
        threads=1,
        workdir=".",
    )
    assert options == expected


# Building pipeline configuration from command line
def args_to_dict(args: argparse.ArgumentParser) -> Dict[str, Any]:
    """
    Parse command line arguments and return a dictionnary ready to be
    dumped into yaml

    Parameters:
        args        ArgumentParser      Parsed arguments from command line

    Return:
                    Dict[str, Any]      A dictionnary containing the parameters
                                        for the pipeline

    Examples:
    >>> example_options = parse_args("/path/to/fasta")
    >>> args_to_dict(example_options)
    {'cold_storage': [' '],
     'design': 'design.tsv',
     'params': {'libType': 'A',
      'salmon_index_extra': '--keepDuplicates --gencode --perfectHash',
      'salmon_quant_extra':
        '--numBootstraps 100 --validateMappings --gcBias --seqBias'},
     'ref': {'fasta': '/path/to/fasta', 'gtf': None},
     'singularity_docker_image': 'docker://continuumio/miniconda3:4.4.10',
     'threads': 1,
     'workdir': '.',
     'workflow': {'aggregate': False, 'fastqc': True, 'multiqc': True}}
    """
    result_dict = {
        "design": os.path.abspath(args.design),
        "config": os.path.abspath(os.path.join(args.workdir, "config.yaml")),
        "workdir": os.path.abspath(args.workdir),
        "threads": args.threads,
        "singularity_docker_image": args.singularity,
        "cold_storage": args.cold_storage,
        "ref": {
            "fasta": os.path.abspath(args.fasta),
            "gtf": (
                os.path.abspath(args.gtf) if args.gtf is not None else None
            ),
        },
        "workflow": {
            "fastqc": not args.no_fastqc,
            "multiqc": not args.no_multiqc,
            "aggregate": args.aggregate,
        },
        "params": {
            "salmon_index_extra": args.salmon_index_extra,
            "salmon_quant_extra": args.salmon_quant_extra,
            "libType": args.libType,
        },
    }
    logging.debug(result_dict)
    return result_dict


def test_args_to_dict() -> None:
    """
    This function simply tests the args_to_dict function with expected output

    Example:
    >>> pytest -v prepare_config.py -k test_args_to_dict
    """
    options = parse_args(
        shlex.split(
            "/path/to/fasta "
            " /path/to/gtf "
            "--design /path/to/design "
            "--workdir /path/to/workdir "
            "--threads 100 "
            "--singularity singularity_image "
            "--cold-storage /path/cold/one /path/cold/two "
            "--no-fastqc "
            "--aggregate "
            "--salmon-index-extra ' --index-arg 1 ' "
            "--salmon-quant-extra ' --quant-arg ok ' "
            "--debug "
        )
    )

    expected = {
        "design": "/path/to/design",
        "config": "/path/to/workdir/config.yaml",
        "workdir": "/path/to/workdir",
        "threads": 100,
        "singularity_docker_image": "singularity_image",
        "cold_storage": ["/path/cold/one", "/path/cold/two"],
        "ref": {"fasta": "/path/to/fasta", "gtf": "/path/to/gtf"},
        "workflow": {"fastqc": False, "multiqc": True, "aggregate": True},
        "params": {
            "salmon_index_extra": " --index-arg 1 ",
            "salmon_quant_extra": " --quant-arg ok ",
            "libType": "A",
        },
    }
    assert args_to_dict(options) == expected


# Yaml formatting
def dict_to_yaml(indict: Dict[str, Any]) -> str:
    """
    This function makes the dictionnary to yaml formatted text

    Parameters:
        indict  Dict[str, Any]  The dictionnary containing the pipeline
                                parameters, extracted from command line

    Return:
                str             The yaml formatted string, directly built
                                from the input dictionnary

    Examples:
    >>> import yaml
    >>> example_dict = {
        "bar": "bar-value",
        "foo": ["foo-list-1", "foo-list-2"]
    }
    >>> dict_to_yaml(example_dict)
    'bar: bar-value\nfoo:\n- foo-list-1\n- foo-list-2\n'
    >>> print(dict_to_yaml(example_dict))
    bar: bar-value
    foo:
    - foo-list-1
    - foo-list-2
    """
    return yaml.dump(indict, default_flow_style=False)


def test_dict_to_yaml() -> None:
    """
    This function tests the dict_to_yaml function with pytest

    Example:
    >>> pytest -v prepare_config.py -k test_dict_to_yaml
    """
    expected = "bar: bar-value\nfoo:\n- foo-list-1\n- foo-list-2\n"
    example_dict = {"bar": "bar-value", "foo": ["foo-list-1", "foo-list-2"]}
    assert dict_to_yaml(example_dict) == expected


# Core of this script
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole configuration sequence

    Parameters:
        args    ArgumentParser      The parsed command line

    Example:
    >>> main(parse_args(shlex.split("/path/to/fasta")))
    """
    # Building pipeline arguments
    logging.debug("Building configuration file:")
    config_params = args_to_dict(args)
    output_path = Path(args.workdir) / "config.yaml"

    # Saving as yaml
    with output_path.open("w") as config_yaml:
        logging.debug(f"Saving results to {str(output_path)}")
        config_yaml.write(dict_to_yaml(config_params))


# Running programm if not imported
if __name__ == "__main__":
    # Parsing command line
    args = parse_args(sys.argv[1:])
    makedirs("logs/prepare")

    # Build logging object and behaviour
    logging.basicConfig(
        filename="logs/prepare/config.log", filemode="w", level=logging.DEBUG
    )

    try:
        logging.debug("Preparing configuration")
        main(args)
    except Exception as e:
        logging.exception("%s", e)
        raise
