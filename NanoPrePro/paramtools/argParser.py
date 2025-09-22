from argparse import ArgumentParser, MetavarTypeHelpFormatter
from NanoPrePro._version import __version__
# TODO: separate each section of the options 
# initiate ArgumentParser
parser = ArgumentParser(
    formatter_class=MetavarTypeHelpFormatter
)

# version
parser.add_argument(
    "--version",
    action="version",
    version=f"%(prog)s {__version__}"
)


# general options
parser.add_argument(
    "--input_fq",
    required=True,
    type=str,
    help="input FASTQ"
)
parser.add_argument(
    "--config",
    type=str,
    help="use the parameters in this config file (JSON)"
    "(can be overriden by command line arguments)"
)
parser.add_argument(
    "--report",
    type=str,
    help="output report file (JSON)",
    default=None
)
parser.add_argument(
    "--processes",
    type=int,
    help="number of processes to use (default: 16)",
    default=16
)
parser.add_argument(
    "--batch_size",
    type=int,
    help="number of records in each batch (default: 1000000)",
    default=1000000
)
parser.add_argument(
    "--seed",
    type=int,
    help="seed for random number generator (default: 42)",
    default=42
)



# optimization options
parser.add_argument(
    "-n",
    type=int,
    help="max number of reads to sample during optimzation (default: 100000)",
    default=100000
)
parser.add_argument(
    "--beta",
    type=float,
    help="the beta parameter for the optimization (default: .1)",
    default=.1
)
parser.add_argument(
    "--min_plen",
    type=float,
    help="the minimal primer length to test (default: .5)",
    default=.5
)


# annotaion options
parser.add_argument(
    "--disable_annot",
    action="store_true",
    help="use this flag to disable annotation"
)
# parser.add_argument(
#     "--skip_lowq",
#     default=0,
#     type=float,
#     help="skip low-quality reads (default: 0)",
# )
# parser.add_argument(
#     "--skip_short",
#     default=0,
#     type=int,
#     help="skip too-short reads (default: 0)"
# )
parser.add_argument(
    "--p5_sense",
    type=str,
    help="5' sense adatper/primer + polyA sequences"
)
parser.add_argument(
    "--p3_sense",
    type=str,
    help="3' sense adatper/primer + polyA sequences"
)
parser.add_argument(
    "--isl5",
    nargs=2,
    type=int,
    help="ideal searching location for 5' adapter/primer sequences (default: optimized)"
)
parser.add_argument(
    "--isl3",
    nargs=2,
    type=int,
    help="ideal searching location for 3' adapter/primer sequences (default: optimized)"
)
parser.add_argument(
    "--pid5",
    type=float,
    help="5' adapter/primer percent identity cutoff (default: optimized)"
)
parser.add_argument(
    "--pid3",
    type=float,
    help="3' adapter/primer percent identity cutoff (default: optimized)"
)
parser.add_argument(
    "--pid_body",
    type=float,
    help="adapter/primer percent identity cutoff (default: optimized)"
)
parser.add_argument(
    "--poly_w",
    type=int,
    help="window size for polyA/T identification",
    default=6
)
parser.add_argument(
    "--poly_k",
    type=int,
    help="number of A/T to be expected in the window",
    default=4
)

# processing (filtering/trimming/orientation) options
parser.add_argument(
    "--trim_adapter",
    action="store_true",
    help="use this flag to trim adatper/primer sequences"
)
parser.add_argument(
    "--trim_poly",
    action="store_true",
    help="use this flag to trim polyA/T sequences"
)
parser.add_argument(
    "--filter_lowq",
    default=7,
    type=float,
    help="filter low-quality reads after all trimming steps (default: 7)"
)
parser.add_argument(
    "--filter_short",
    default=0,
    type=int,
    help="filter too short reads after all trimming steps (default: 0)"
)
parser.add_argument(
    "--orientation",
    default=1,
    type=int,
    help="re-orient reads (0: generic (default), 1: sense, -1: antisense)"
)

# output options
parser.add_argument(
    "--output_fusion",
    type=str,
    help="output fusion/chimeric reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--output_truncated",
    type=str,
    help="output truncated/non-full-length reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--output_full_length",
    type=str,
    help="output full-length reads to this file (use '-' for stdout)",
    default="-"
)
parser.add_argument(
    "--suffix_filtered",
    type=str,
    help="output filtered reads with the suffix"
)
