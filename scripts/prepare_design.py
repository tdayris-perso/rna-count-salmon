#!/usr/bin/python3.7
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
"""

import argparse                      # Parse command line
import logging                       # Traces and loggings
import logging.handlers              # Logging behaviour
import os                            # OS related activities
import pandas as pd                  # Parse TSV files
from pathlib import Path             # Paths related methods
import pytest                        # Unit testing
import shlex                         # Lexical analysis
import sys                           # System related methods
from typing import Dict, Generator, List   # Type hints

logger = logging.getLogger(
    os.path.splitext(os.path.basename(sys.argv[0]))[0]
)


# Building custom class for help formatter
class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """
    pass


# Processing functions
# Looking for fastq files
def search_fq(fq_dir: Path,
              recursive: bool = False) -> Generator[str, str, None]:
    """
    Iterate over a directory and search for fastq files
    """
    for path in fq_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_fq(path, recursive)
            else:
                continue

        if path.name.endswith(("fq", "fq.gz", "fastq", "fastq.gz")):
            yield path


# Testing search_fq
def test_search_fq():
    """
    This function tests the ability of the function "search_fq" to find the
    fastq files in the given directory
    """
    path = Path("../tests/reads/")
    expected = list(
        Path("../tests/reads/{}_R{}.fastq".format(sample, stream))
        for sample in ["A", "B"]
        for stream in [1, 2]
    )
    assert sorted(list(search_fq(path))) == sorted(expected)


# Handling logging options
# No tests for this function
def setup_logging(args: argparse.ArgumentParser) -> None:
    """
    Configure logging behaviour
    """
    root = logging.getLogger("")
    root.setLevel(logging.WARNING)
    logger.setLevel(args.debug and logging.DEBUG or logging.INFO)
    if not args.quiet:
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(
            "%(levelname)s [%(name)s]: %(message)s"
        ))
        root.addHandler(ch)


# Turning the FQ list into a dictionnary
def classify_fq(fq_files: List[Path], paired: bool = True) -> Dict[str, Path]:
    """
    Return a dictionnary with identified fastq files (paried or not)
    """
    fq_dict = {}
    if paired is not True:
        # Case single fastq per sample
        logger.debug("Sorting fastq files as single-ended")
        for fq in fq_files:
            fq_dict[fq.name] = {
                "Sample_id": fq.stem,
                "Upstream_file": fq.absolute()
            }
    else:
        # Case pairs of fastq are used
        logger.debug("Sorting fastq files as pair-ended")
        for fq1, fq2 in zip(fq_files[0::2], fq_files[1::2]):
            fq_dict[fq1.name] = {
                "Sample_id": fq1.stem,
                "Upstream_file": fq1.absolute(),
                "Downstream_file": fq2.absolute()
            }
    print(fq_dict)
    return fq_dict


# Testing the classification function
def test_classify_fq():
    """
    This function takes input from the pytest decorator
    to test the classify_fq function
    """
    prefix = Path(__file__).parent.parent
    fq_list = sorted(list(search_fq(prefix / "tests" / "reads")))
    expected = {
        'A_R1.fastq': {
            'Sample_id': 'A_R1',
            'Upstream_file': prefix / "tests" / "reads" / "A_R1.fastq",
            'Downstream_file': prefix / "tests" / "reads" / 'A_R2.fastq'
        },
        'B_R1.fastq': {
            'Sample_id': 'B_R1',
            'Upstream_file': prefix / "tests" / "reads" / 'B_R1.fastq',
            'Downstream_file': prefix / "tests" / "reads" / 'B_R2.fastq'
        }
    }
    logger = logging.getLogger(
        os.path.splitext(os.path.basename(sys.argv[0]))[0]
    )
    assert classify_fq(fq_list) == expected


# Parsing command line arguments
# This function won't be tested
def parse_args(args: str = sys.argv[1:]) -> argparse.ArgumentParser:
    """
    Simply build a command line parser object
    """
    # Defining command line options
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter
    )

    # Required arguments
    main_parser.add_argument(
        "path",
        help="Path to the directory containing fastq files",
        type=str
    )

    # Optional arguments
    main_parser.add_argument(
        "-s", "--single",
        help="The samples are single ended rnaseq reads, not pair ended",
        action="store_true"
    )

    main_parser.add_argument(
        "-r", "--recursive",
        help="Recursively search in sub-directories for fastq files",
        action="store_true"
    )

    main_parser.add_argument(
        "-o", "--output",
        help="Path to output file (default: %(default)s)",
        type=str,
        default="design.tsv"
    )

    # Logging options
    log = main_parser.add_mutually_exclusive_group()
    log.add_argument(
        "-d", "--debug",
        help="Set logging in debug mode",
        default=False,
        action='store_true'
    )

    log.add_argument(
        "-q", "--quiet",
        help="Turn off logging behaviour",
        default=False,
        action='store_true'
    )

    # Parsing command lines
    return main_parser.parse_args()


# Main function, the core of this programm
def main(args: argparse.ArgumentParser) -> None:
    """
    This function performs the whole preparation sequence
    """
    # Searching for fastq files and sorting them alphabetically
    fq_files = sorted(list(search_fq(Path(args.path), args.recursive)))
    logger.debug("Head of alphabeticaly sorted list of fastq files:")
    logger.debug([str(i) for i in fq_files[0:5]])

    # Building a dictionnary of fastq (pairs?) and identifiers
    fq_dict = classify_fq(fq_files)

    # Using Pandas to handle TSV output (yes pretty harsh I know)
    data = pd.DataFrame(fq_dict).T
    logger.debug("\n{}".format(data.head()))
    logger.debug("Saving results to {}".format(args.output))
    data.to_csv(args.output, sep="\t", index=False)


# Running programm if not imported
if __name__ == '__main__':
    args = parse_args()
    logger = logging.getLogger(
        os.path.splitext(os.path.basename(sys.argv[0]))[0]
    )
    setup_logging(args)

    try:
        logger.debug("Preparing design")
        main(args)
    except Exception as e:
        logger.exception("%s", e)
        sys.exit(1)
    sys.exit(0)
