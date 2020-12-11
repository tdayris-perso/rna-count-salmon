#!/usr/bin/python3.8
# -*- coding: utf-8 -*-


"""
This is the CLI launcher of the rna-count-salmon pipeline. This pipeline is
powered by Snakemake and uses Salmon and FastQC on raw fastq files, then
MultiQC on the results to aggregate quality reports.

You will have quantified reads over both transcriptomic and genomic regions,
alongside with complete quality reports. Please be aware that a Snakemake report
is available for more details on both results and methods.

If you have any question, please refer to the wiki at:
https://github.com/tdayris-perso/rna-count-salmon/wiki
Or open an issue at:
https://github.com/tdayris-perso/rna-count-salmon/issues

Citations:
https://github.com/tdayris-perso/rna-count-salmon/wiki/Pipeline_Content
"""

import argparse
import pytest
import os
import logging
import sys

from pathlib import Path
from snakemake.utils import makedirs
from snakemake.shell import shell
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
        formatter_class=common_script_rna_count_salmon.CustomFormatter,
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
        "--snakemake-args",
        help="Snakemake arguments. If you use this wrapper instead "
             "of real snakemake call, then please put all your"
             " arguments in simple quotes. (default: %(default)s)",
        type=str,
        default=""
    )
    snake.add_argument(
        "--no-profile",
        help="Do not activate snakemake profile. This means you have to "
             "define several environment variables. See documentation for "
             "more information.",
        action="store_true",
        default=False
    )
    snake.add_argument(
        "--cache",
        help="Use the shared cache that allows user not to re-run "
             "indexation steps. Not using cache will slow down your analysis "
             "and use disk space for duplicated data.",
        action="store_true",
        default=False
    )
    snake.set_defaults(func=snakemake_run)

    report_parser = subparsers.add_parser(
        "report",
        add_help=True
    )
    report_parser.add_argument(
        "--snakemake-args",
        help="Snakemake arguments. If you use this wrapper instead "
             "of real snakemake call, then please put all your"
             " arguments in simple quotes. (default: %(default)s)",
        type=str,
        default=""
    )
    report_parser.add_argument(
        "--no-profile",
        help="Do not activate snakemake profile. This means you have to "
             "define several environment variables. See documentation for "
             "more information.",
        action="store_true",
        default=False
    )
    report_parser.add_argument(
        "--cache",
        help="Use the shared cache that allows user not to re-run "
             "indexation steps. Not using cache will slow down your analysis "
             "and use disk space for duplicated data.",
        action="store_true",
        default=False
    )

    igr = subparsers.add_parser(
        "flamingo",
        add_help=True
    )
    igr.add_argument(
        "--fastqdir",
        help="Path to fastq files directory (default: %(default)s)",
        type=str,
        metavar="FQ-DIR",
        default=os.getcwd()
    )
    igr.set_defaults(func=igr_run)

    report_parser.set_defaults(func=report)
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


def snakemake_command(opt: str = "",
                      use_profile: bool = True,
                      make_report: bool = False,
                      use_cache: bool = True) -> str:
    """
    Build snakemake command line
    """
    return [
        "snakemake",
        f"-s {os.getenv('SNAKEFILE')}",
        f"--profile {os.getenv('PROFILE')}" if use_profile is True else "",
        "--report Quantification_Report.html" if make_report is True else "",
        "--cache salmon_index tr2gene" if use_cache is True else "",
        opt
    ]


@pytest.mark.parametrize(
    "opt, prof, rep, cache, expected", [
        ("--configfile config.yaml", True, True, True, [
            "snakemake",
            f"-s {os.getenv('SNAKEFILE')}",
            f"--profile {os.getenv('PROFILE')}",
            "--report Quantification_Report.html",
            "--cache salmon_index tr2gene",
            "--configfile config.yaml"
        ]),
        ("", False, False, False, [
            "snakemake", f"-s {os.getenv('SNAKEFILE')}", "", "", "", ""
        ]),
    ]
)
def test_snakemake_command(opt, prof, rep, cache, expected) -> None:
    """
    This function tests the above snakemake_command
    """
    assert snakemake_command(opt, prof, rep, cache) == expected


def check_env() -> bool:
    """
    Verify environment variables required for this pipeline
    """
    expected = [
        os.getenv('SNAKEFILE'),
        os.getenv('PROFILE'),
        os.getenv("RNA_COUNT_LAUNCHER"),
        os.getenv("FASTA"),
        os.getenv("GTF")
    ]
    test_if_none = all(var is not None for var in expected)

    try:
        test_if_exists = all(os.path.exists(path) for path in expected)
    except TypeError:
        print(f"Environment does not fit expectations: {expected}")

    return test_if_none and test_if_exists



def run_cmd(*cmd_line) -> None:
    """
    Run a provided command line
    """
    cmd_line = " ".join(cmd_line)
    if check_env():
        print(cmd_line)
        shell(cmd_line)
    else:
        print("Environment was not suitable for this pipeline to run.")

def snakemake_run(cmd_line_args) -> None:
    """
    Call snakemake itself
    """
    command = snakemake_command(
        opt=cmd_line_args.snakemake_args,
        use_profile=not cmd_line_args.no_profile,
        make_report=False,
        use_cache=cmd_line_args.cache
    )

    run_cmd(*command)


def report(cmd_line_args) -> None:
    """
    Call snakemake itself
    """
    command = snakemake_command(
        opt=cmd_line_args.snakemake_args,
        use_profile=not cmd_line_args.no_profile,
        make_report=True,
        use_cache=cmd_line_args.cache
    )

    run_cmd(*command)


def igr_run(cmd_line_args) -> None:
    """
    Call this pipeline whole pipeline with default arguments
    """
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        config_cmd = [
            "python3",
            os.getenv('RNA_COUNT_LAUNCHER'),
            "config",
            os.getenv('FASTA'),
            os.getenv('GTF'),
            "--threads 20",
            "--aggregate",
            "--debug",
            "--cold-storage /mnt/isilon /mnt/archivage",
        ]
        run_cmd(*config_cmd)

    design_path = "design.tsv"
    if not os.path.exists(design_path):
        design_cmd = [
            "python3",
            os.getenv('RNA_COUNT_LAUNCHER'),
            "design",
            cmd_line_args.fastqdir,
            "--recursive",
            "--debug"
        ]
        run_cmd(*design_cmd)

    snakemake_cmd = ["python3", os.getenv("RNA_COUNT_LAUNCHER"), "snakemake"]
    run_cmd(*snakemake_cmd)

    report_cmd = ["python3", os.getenv("RNA_COUNT_LAUNCHER"), "report"]
    run_cmd(*report_cmd)



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
