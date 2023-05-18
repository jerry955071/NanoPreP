from argparse import ArgumentParser

# initiate ArgumentParser
parser = ArgumentParser(description='Arguments availible to nanoprep-optimize')

# positional arguments
parser.add_argument(
    "input_fq",
    metavar="input.fq",
    type=str,
    help="input FASTQ"
)

# general options
parser.add_argument(
    "-n",
    metavar="int",
    type=int,
    help="sample `n` reads from `input_fq` (default: 10000)"
)
parser.add_argument(
    "--config",
    metavar="PATH",
    type=str,
    help="use the parameters in this config file (JSON) "
    "(can be overriden by command line arguments)"
)

# output options
parser.add_argument(
    "--output",
    metavar="PATH",
    type=str,
    help="output html"
)
# parser.add_argument(
#     "--output_config",
#     metavar="PATH",
#     default="config.json",
#     type=str,
#     help="output config.json (default: config.json)"
# )

# annotation options
parser.add_argument(
    "--skip_lowq",
    metavar="float",
    type=float,
    help="skip low-quality reads (default: -1)"
)
parser.add_argument(
    "--skip_short",
    metavar="int",
    type=int,
    help="skip too-short reads (default: -1)"
)
parser.add_argument(
    "--p5_sense",
    metavar="str",
    type=str,
    help="5' sense adatper/primer sequences"
)
parser.add_argument(
    "--p3_sense",
    metavar="str",
    type=str,
    help="3' sense adatper/primer sequences"
)
parser.add_argument(
    "--isl5",
    metavar="int",
    nargs=2,
    type=int,
    help="ideal searching location for 5' adapter/primer sequences "
    "(e.g. 1 130)"
)
parser.add_argument(
    "--isl3",
    metavar="int",
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
    "skip_lowq": 7,
    "skip_short": 190,
    "p5_sense": "GTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTTT",
    "p3_sense": "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAG",
    "isl5": [0, 130],
    "isl3": [-60, -1]
}