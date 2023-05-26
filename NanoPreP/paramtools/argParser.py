from argparse import ArgumentParser, MetavarTypeHelpFormatter
from NanoPreP._version import __version__
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
    "--mode",
    type=str,
    help="use presets standard/annotate/report "
    "(can be overriden by `config` and command line arguments) "
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
    help="output report file (JSON)"
)

# optimization options
parser.add_argument(
    "-n",
    type=int,
    help="max number of reads to sample during optimzation (default: 100000)"
)
parser.add_argument(
    "--precision",
    type=float,
    help="precision cutoff while optimizing pid cutoff"
)

# annotaion options
parser.add_argument(
    "--disable_annot",
    action="store_true",
    help="use this flag to disable annotation"
)
parser.add_argument(
    "--skip_lowq",
    default=-1,
    type=float,
    help="skip low-quality reads (default: -1)"
)
parser.add_argument(
    "--skip_short",
    default=-1,
    type=int,
    help="skip too-short reads (default: -1)"
)
parser.add_argument(
    "--p5_sense",
    type=str,
    help="5' sense adatper/primer/polymer sequences"
)
parser.add_argument(
    "--p3_sense",
    type=str,
    help="3' sense adatper/primer/polymer sequences"
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
parser.add_argument(
    "--pid_isl",
    type=float,
    help="adapter/primer percent identity cutoff (in ISL)"
)
parser.add_argument(
    "--pid_body",
    type=float,
    help="adapter/primer percent identity cutoff (on read body)"
)
parser.add_argument(
    "--poly_w",
    type=int,
    help="window size for homopolymer identification"
)
parser.add_argument(
    "--poly_k",
    type=int,
    help="number of monomers to be expected in the window"
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
    help="use this flag to trim homopolymers"
)
parser.add_argument(
    "--filter_lowq",
    default=-1,
    type=float,
    help="filter low-quality reads after all trimming steps (default: -1)"
)
parser.add_argument(
    "--filter_short",
    default=-1,
    type=int,
    help="filter too short reads after all trimming steps (default: -1)"
)
parser.add_argument(
    "--orientation",
    default=0,
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
    help="output full-length reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--suffix_filtered",
    type=str,
    help="output filtered reads with the suffix"
)
