# -*- coding: UTF-8 -*-

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLv3"


import sys

if sys.version_info < (3, 5):
    raise ValueError("Python 3.5 is required.\n", file=sys.stderr)


try:
    import setuptools
except ImportError:
    raise ImportError("Please install setuptools.", file=sys.stderr)

requirements = [
    "conda",
    "datrie ==0.8.2",
    "jinja2 ==2.11.2",
    "flask ==1.1.2",
    "pandas ==1.0.3",
    "snakemake ==5.19.3",
    "pytest ==5.4.2",
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rna-count-salmon",
    version="v1.9",
    author="Thibault Dayris",
    author_email="thibault.dayris@gustaveroussy.fr",
    description="rna-count-salmon is a pipeline written with "
                " Snakemake that quantifies your RNA-Seq reads "
                "with Salmon, controls input quality with FastQC, "
                "and gathers qualities with MultiQC.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    zip_safe=False,
    license="GPLv3",
    url="https://github.com/tdayris-perso/rna-count-salmon/wiki",
    install_requires=requirements,
    test_requires=requirements,
    setup_requires=requirements,
    packages=setuptools.find_packages(),
    packages_data={
        "snakemake_rules": setuptools.find_packages("rules", include="*.smk|py"),
        "snakemake_scripts": setuptools.find_packages("scripts", include="*.py"),
        "snakemake": "Snakefile"
    },
    scripts=["rna-count-salmon.py"],
    python_requires=">=3.7",
    classifiers=[
        "Environment :: Console",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 5 - Production/Stable"
    ]
)
