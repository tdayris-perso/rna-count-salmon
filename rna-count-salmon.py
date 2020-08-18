#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This is the CLI launcher of the rna-count-salmon pipeline. This pipeline is
powered by Snakemake and uses Salmon and FastQC on raw fastq files, then
MultiQC on the results to aggregate quality reports.

You will have quantified reads over both transcriptomic and genomic regions,
alongside with complete quality reports.

Please be aware that a Snakemake report is available for more details on
both results and methods.

If you have any question, please refer to the wiki at:
https://github.com/tdayris-perso/rna-count-salmon/wiki

Or open an issue at:
https://github.com/tdayris-perso/rna-count-salmon/issues

Citations:
https://github.com/tdayris-perso/rna-count-salmon/wiki/Pipeline_Content
"""

import argparse
import os
import logging
import sys

from pathlib import Path
from snakemake.utils import makedirs
from snakemake.shell import shell
import snakemake
from typing import Any

try:
    from scripts import prepare_config, prepare_design, common_script_rna_count_salmon
except ImportError:
    scripts_path = Path(os.path.realpath(__file__)).parent / "scripts"
    sys.path.append(str(scripts_path))
    import prepare_config, prepare_design, common_script_rna_count_salmon
except ModuleNotFoundError:
    scripts_path = Path(os.path.realpath(__file__)).parent / "scripts"
    sys.path.append(str(scripts_path))
    import prepare_config, prepare_design, common_script_rna_count_salmon



def parser() -> argparse.ArgumentParser:
    """
    Build command line parser object
    """
    main_parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
    )

    subparsers = main_parser.add_subparsers()
    config = subparsers.add_parser(
        "config",
        parents=[prepare_config.parser()],
        add_help=False
    )
    config.set_defaults(func=prepare_config.main)

    design = subparsers.add_parser(
        "design",
        parents=[prepare_design.parser()],
        add_help=False
    )
    design.set_defaults(func=prepare_design.main)

    snake = subparsers.add_parser(
        "snakemake",
        add_help=True
    )
    snake.add_argument(
        "ARGS",
        help="Snakemake arguments. If you use this wrapper instead "
             "of real snakemake call, then please put all your"
             " arguments in simple quotes.",
        type=str,
        nargs="+",
        default="--help"
    )
    snake.set_defaults(func=snakefile)
    return main_parser


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


def snakefile(args):
    """
    Call snakemake itself
    """
    snakefile_path = Path(os.path.realpath(__file__))
    snakefile_path = snakefile_path.parent / "Snakefile"

    command = f"snakemake -s {snakefile_path} {' '.join(args.ARGS)}"

    shell(command)


if __name__ == '__main__':
    # Parsing command line arguments
    args = parse_args(sys.argv[1:])
    makedirs("logs/prepare")

    # Build logging object and behaviour
    logging.basicConfig(
        filename="logs/prepare/design.log", filemode="w", level=logging.DEBUG
    )

    try:
        args.func(args)
    except Exception as e:
        logging.exception("%s", e)
        raise
