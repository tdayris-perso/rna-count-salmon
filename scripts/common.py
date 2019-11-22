#!/usr/bin/python3.7
# -*- coding: utf-8 -*-


"""
This script contains functions that are to be called by any other scripts in
this pipeline.
"""

import argparse    # Argument parsing
import logging     # Logging behaviour


# Building custom class for help formatter
class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    """
    This class is used only to allow line breaks in the documentation,
    without breaking the classic argument formatting.
    """
    pass


# Handling logging options
# No tests for this function
def setup_logging(logger: logging.Logger,
                  args: argparse.ArgumentParser = None) -> None:
    """
    Configure logging behaviour
    """
    root = logging.getLogger("")
    root.setLevel(logging.WARNING)
    logger.setLevel(args.debug and logging.DEBUG or logging.INFO)
    if (args is None) or (args.quiet is False):
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(
            "%(levelname)s [%(name)s]: %(message)s"
        ))
        root.addHandler(ch)