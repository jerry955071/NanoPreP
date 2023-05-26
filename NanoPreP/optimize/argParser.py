"""Argument parser for NanoPreP optimize"""
from argparse import ArgumentParser, MetavarTypeHelpFormatter

# initiate ArgumentParser
parser = ArgumentParser(
    formatter_class=MetavarTypeHelpFormatter
)

# general options
parser.add_argument(
    "--input_fq",
    required=True,
    type=str,
    help="input FASTQ"
)
parser.add_argument(
    "-n",
    type=int,
    help="max number of reads to sample during optimzation (default: 100000)"
)
parser.add_argument(
    "--config",
    type=str,
    help="use the parameters in this config file (JSON) "
    "(can be overriden by command line arguments)"
)

# output options
parser.add_argument(
    "--output",
    type=str,
    help="output html"
)

# annotation options
parser.add_argument(
    "--skip_lowq",
    type=float,
    help="skip low-quality reads (default: -1)"
)
parser.add_argument(
    "--skip_short",
    type=int,
    help="skip too-short reads (default: -1)"
)
parser.add_argument(
    "--p5_sense",
    type=str,
    help="5' sense adatper/primer sequences"
)
parser.add_argument(
    "--p3_sense",
    type=str,
    help="3' sense adatper/primer sequences"
)
parser.add_argument(
    "--isl5",
    nargs=2,
    type=int,
    help="ideal searching location for 5' adapter/primer sequences "
    "(e.g. 0 130)"
)
parser.add_argument(
    "--isl3",
    nargs=2,
    type=int,
    help="ideal searching location for 3' adapter/primer sequences "
    "(e.g. -60 -1)"
)

# Default parameters
Defaults = {
    "input_fq": None,
    "n": 100000,
    "config": None,
    "output": "output.html",
    "skip_lowq": -1,
    "skip_short": -1,
    "p5_sense": "TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG",
    "p3_sense": "GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC",
    "isl5": [0, 130],
    "isl3": [-60, -1],
    "precision": 0.99
}