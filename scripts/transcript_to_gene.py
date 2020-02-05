#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script extracts a tsv from a GTF. This tsv contains gene_id,
transcript_id, gene_name. This is used by the aggregation step to
annotate trnascript identifiers with human readable gene names.

Input: GTF
Output:

- Column 1: Gene identifier
- Column 2: Transcript identifier
- Column 3: Gene name
- Column 4: Chromorome
- Column 5: Transcript start
- Column 6: Transcript end
- Column 7: Transcript strand
"""

import pytest

from pathlib import Path
from typing import List, Optional

# Pytest variables
short_comment = "#test"
long_comment = '#21\tok\ttranscript\t145\t136\t.\t+\t0\tgene_id "ge_ok"'
gtf_example1 = (
    '21\tok\ttranscript\t145\t136\t.\t+\t0'
    '\tgene_id "ge_ok"; transcript_id "tr_ok"'
)
gtf_example2 = (
    '21\tok\ttranscript\t145\t136\t.\t+\t0\t'
    'gene_id "ge_ok"; transcript_id "tr_ok"; gene_name "name_ok"'
)


def parse_gtf(line: str) -> Optional[List[str]]:
    """
    This function take a gtf line as input and returns either None if the line
    does not contain any transcript information, or the list described in the
    general docstring above.
    """
    if line.startswith("#"):
        # Only comments start with '#'
        # they shall be ignored
        return None

    # If the line does not describe a transcript, the we do not care
    # about it, thus it would not contain all the required information
    chomp = line[:-1].split("\t")
    if chomp[2] != "transcript":
        return None

    start = chomp[3]
    end = chomp[4]
    chrom = chomp[0]
    strand = chomp[6]

    # Parse the last column
    chomp = {
        attr.split('"')[0].strip(): attr.split('"')[1].strip()
        for attr in chomp[8].split(";")
        if attr != '' and '"' in attr
    }

    # Some genes have an ID but no name ...
    try:
        result = [
            chomp["gene_id"],
            chomp["transcript_id"],
            chomp["gene_name"],
            chrom,
            start,
            end,
            strand
        ]
    except KeyError:
        result = [
            chomp["gene_id"],
            chomp["transcript_id"],
            chomp["gene_id"],
            chrom,
            start,
            end,
            strand
        ]

    return result


@pytest.mark.parametrize(
    "line, expected", [
        (short_comment, None),
        (long_comment, None),
        (gtf_example1, ["ge_ok", "tr_ok", "ge_ok", "21", "145", "136", "+"]),
        (gtf_example2, ["ge_ok", "tr_ok", "name_ok", "21", "145", "136", "+"])
    ]
)
def test_parse_gtf(line: str, expected: Optional[List[str]]) -> None:
    """
    This function simply tests the gtf parsing with multiple arguments
    """
    assert parse_gtf(line) == expected


# Main programm ran only if ran by main (snakemake//pytest protection)
if __name__ == '__main__':
    gtf_path = Path(snakemake.input["gtf"])
    output_path = Path(snakemake.output["tsv"])

    with gtf_path.open("r") as gtf, output_path.open("w") as tsv:
        print(
            "\t".join(
                ["gene_id", "transcript_id", "gene_name",
                 "chrom", "start", "end", "strand"]
            ),
            file=tsv
        )

        for line in gtf:
            line = parse_gtf(line)
            if line is not None:
                print("\t".join(line), file=tsv)
