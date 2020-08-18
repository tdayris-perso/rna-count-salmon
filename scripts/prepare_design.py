#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
This script aims to prepare the list of files to be processed
by the rna-count-salmon pipeline

It iterates over a given directory, lists all fastq files. As a pair of
fastq files usually have names that follows each other in the alphabetical
order, this script sorts the fastq files names and, by default, creates
pairs of fastq files that way.

Finally, it writes these pairs, using the longest common substring as
identifier. The written file is a TSV file.

You can test this script with:
pytest -v ./prepare_design.py

Usage example:
# Single ended reads example:
python3.7 ./prepare_design.py tests/reads --single

# Paired-end libary example:
python3.7 ./prepare_design.py tests/reads

# Search in sub-directories:
python3.7 ./prepare_design.py tests --recursive
"""

import argparse  # Parse command line
import logging  # Traces and loggings
import os  # OS related activities
import pandas as pd  # Parse TSV files
import pytest  # Unit testing
import shlex  # Lexical analysis
import sys  # System related methods

from pathlib import Path  # Paths related methods
from snakemake.utils import makedirs  # Easily build directories
from typing import Dict, Generator, List, Any  # Type hints

try:
    from scripts.common_script_rna_count_salmon import *
except ModuleNotFoundError:
    from common_script_rna_count_salmon import *


# Processing functions
# Looking for fastq files
def search_fq(
    fq_dir: Path, recursive: bool = False
) -> Generator[str, str, None]:
    """
    Iterate over a directory and search for fastq files

    Parameters:
        fq_dir      Path        Path to the fastq directory in which to search
        recursive   bool        A boolean, weather to search recursively in
                                sub-directories (True) or not (False)

    Return:
                    Generator[str, str, None]       A Generator of paths

    Example:
    >>> search_fq(Path("tests/reads/"))
    <generator object search_fq at 0xXXXXXXXXXXXX>

    >>> list(search_fq(Path("tests/", True)))
    [PosixPath('tests/reads/A_R2.fastq'),
     PosixPath('tests/reads/B_R2.fastq'),
     PosixPath('tests/reads/A_R1.fastq'),
     PosixPath('tests/reads/B_R1.fastq')]
    """
    for path in fq_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_fq(path, recursive)
            else:
                continue

        if path.name.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
            yield path


# Testing search_fq
def test_search_fq():
    """
    This function tests the ability of the function "search_fq" to find the
    fastq files in the given directory

    Example:
    pytest -v prepare_design.py -k test_search_fq
    """
    path = Path("tests/reads/")

    expected = list(
        path / "{}_R{}.fastq".format(sample, stream)
        for sample in ["A", "B"]
        for stream in [1, 2]
    )
    assert sorted(list(search_fq(path))) == sorted(expected)


# Turning the FQ list into a dictionnary
def classify_fq(fq_files: List[Path], paired: bool = True) -> Dict[str, Path]:
    """
    Return a dictionnary with identified fastq files (paried or not)

    Parameters:
        fq_files    List[Path]      A list of paths to iterate over
        paired      bool            A boolean, weather the dataset is
                                    pair-ended (True) or single-ended (False)

    Return:
                    Dict[str, Path] A dictionnary: for each Sample ID, the ID
                                    is repeated alongside with the upstream
                                    /downstream fastq files.

    Example:
    # Paired-end single sample
    >>> classify_fq([Path("file1.R1.fq"), Path("file1.R2.fq")], True)
    {'file1.R1.fq': {'Downstream_file': PosixPath("/path/to/file1.R1.fq"),
     'Sample_id': 'file1.R1',
     'Upstream_file': PosixPath('/path/to/file1.R2.fq')}

    # Single-ended single sample
    >>> classify_fq([Path("file1.fq")], False)
    {'file1.fq': {'Sample_id': 'file1',
     'Upstream_file': PosixPath('/path/to/file1.fq')}}
    """
    fq_dict = {}
    if paired is not True:
        # Case single fastq per sample
        logging.debug("Sorting fastq files as single-ended")
        for fq in fq_files:
            fq_dict[fq.name] = {
                "Sample_id": fq.stem,
                "Upstream_file": fq.absolute(),
            }
    else:
        # Case pairs of fastq are used
        logging.debug("Sorting fastq files as pair-ended")
        for fq1, fq2 in zip(fq_files[0::2], fq_files[1::2]):
            fq_dict[fq1.name] = {
                "Sample_id": fq1.stem,
                "Upstream_file": fq1.absolute(),
                "Downstream_file": fq2.absolute(),
            }
    logging.debug(fq_dict)
    return fq_dict


def test_classify_fq():
    """
    This function takes input from the pytest decorator
    to test the classify_fq function

    Example:
    pytest -v ./prepare_design.py -k test_classify_fq
    """
    prefix = Path(__file__).parent.parent
    fq_list = sorted(list(search_fq(prefix / "tests" / "reads")))
    expected = {
        "A_R1.fastq": {
            "Sample_id": "A_R1",
            "Upstream_file": prefix / "tests" / "reads" / "A_R1.fastq",
            "Downstream_file": prefix / "tests" / "reads" / "A_R2.fastq",
        },
        "B_R1.fastq": {
            "Sample_id": "B_R1",
            "Upstream_file": prefix / "tests" / "reads" / "B_R1.fastq",
            "Downstream_file": prefix / "tests" / "reads" / "B_R2.fastq",
        },
    }

    assert classify_fq(fq_list) == expected


# Parsing command line arguments
# This function won't be tested
def parser() -> argparse.ArgumentParser:
    """
    Build a command line parser object

    Parameters:
        args    Any                 Command line arguments

    Return:
                ArgumentParser      Parsed command line object
    """
    # Defining command line options
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter,
        epilog="This script does not perform any magic. Check the result.",
    )

    # Required arguments
    main_parser.add_argument(
        "path",
        help="Path to the directory containing fastq files",
        type=str
    )

    # Optional arguments
    main_parser.add_argument(
        "-s",
        "--single",
        help="The samples are single ended rnaseq reads, not pair ended",
        action="store_true",
    )

    main_parser.add_argument(
        "-r",
        "--recursive",
        help="Recursively search in sub-directories for fastq files",
        action="store_true",
    )

    main_parser.add_argument(
        "-o",
        "--output",
        help="Path to output file (default: %(default)s)",
        type=str,
        default="design.tsv",
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

    # Parsing command lines
    return main_parser


def parse_args(args: Any = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    Parse command line arguments

    Parameters:
        args    Any                 Command line arguments

    Return:
                ArgumentParser      Parsed command line object

    Example:
    >>> parse_args(shlex.split("/path/to/fasta --single"))
    Namespace(debug=False, output='design.tsv', path='/path/to/fasta',
    quiet=False, recursive=False, single=True)
    """
    # Parsing command lines
    return parser().parse_args(args)


def test_parse_args() -> None:
    """
    This function tests the command line parsing

    Example:
    >>> pytest -v prepare_config.py -k test_parse_args
    """
    options = parse_args(shlex.split("/path/to/fastq/dir/"))
    expected = argparse.Namespace(
        debug=False,
        output="design.tsv",
        path="/path/to/fastq/dir/",
        quiet=False,
        recursive=False,
        single=False,
    )
    assert options == expected


# Main function, the core of this script
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole preparation sequence

    Parameters:
        args    ArgumentParser      The parsed command line

    Example:
    >>> main(parse_args(shlex.split("/path/to/fasta/dir/")))
    """
    # Searching for fastq files and sorting them alphabetically
    fq_files = sorted(list(search_fq(Path(args.path), args.recursive)))
    logging.debug("Head of alphabeticaly sorted list of fastq files:")
    logging.debug([str(i) for i in fq_files[0:5]])

    # Building a dictionnary of fastq (pairs?) and identifiers
    fq_dict = classify_fq(fq_files)

    # Using Pandas to handle TSV output (yes pretty harsh I know)
    data = pd.DataFrame(fq_dict).T
    logging.debug("\n{}".format(data.head()))
    logging.debug("Saving results to {}".format(args.output))
    data.to_csv(args.output, sep="\t", index=False)


# Running programm if not imported
if __name__ == "__main__":
    # Parsing command line
    args = parse_args(sys.argv[1:])

    makedirs("logs/prepare")

    # Build logging object and behaviour
    logging.basicConfig(
        filename="logs/prepare/design.log", filemode="w", level=logging.DEBUG
    )

    try:
        logging.debug("Preparing design")
        main(args)
    except Exception as e:
        logging.exception("%s", e)
        raise
