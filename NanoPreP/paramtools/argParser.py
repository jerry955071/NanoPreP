from argparse import ArgumentParser
# TODO: separate each section of the options 
# initiate ArgumentParser
parser = ArgumentParser(description='Arguments availible to NanoPreP')

# positional arguments
parser.add_argument(
    "input_fq",
    metavar="input.fq",
    type=str,
    help="input FASTQ"
)

# general options
parser.add_argument(
    "--mode",
    metavar="[strandard|annotate|report]",
    type=str,
    help="use parameter presets "
    "(can be overriden by `config` and command line arguments) "
)
parser.add_argument(
    "--config",
    metavar="PATH",
    type=str,
    help="use the parameters in this config file (JSON)"
    "(can be overriden by command line arguments)"
)
parser.add_argument(
    "--report",
    metavar="PATH",
    type=str,
    help="output report file (JSON)"
)

# annotaion options
parser.add_argument(
    "--disable_annot",
    action="store_true",
    help="use this flag to disable annotation"
)
parser.add_argument(
    "--skip_lowq",
    metavar="float",
    default=-1,
    type=float,
    help="skip low-quality reads (default: -1)"
)
parser.add_argument(
    "--skip_short",
    metavar="int",
    default=-1,
    type=int,
    help="skip too-short reads (default: -1)"
)
parser.add_argument(
    "--p5_sense",
    metavar="str",
    type=str,
    help="5' sense adatper/primer/polymer sequences"
)
parser.add_argument(
    "--p3_sense",
    metavar="str",
    type=str,
    help="3' sense adatper/primer/polymer sequences"
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
parser.add_argument(
    "--pid_isl",
    metavar="float",
    type=float,
    help="adapter/primer percent identity cutoff (in ISL)"
)
parser.add_argument(
    "--pid_body",
    metavar="float",
    type=float,
    help="adapter/primer percent identity cutoff (on read body)"
)
parser.add_argument(
    "--poly_w",
    metavar="int",
    type=int,
    help="window size for homopolymer identification"
)
parser.add_argument(
    "--poly_k",
    metavar="int",
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
    metavar="float",
    default=-1,
    type=float,
    help="filter low-quality reads after all trimming steps (default: -1)"
)
parser.add_argument(
    "--filter_short",
    metavar="int",
    default=-1,
    type=int,
    help="filter too short reads after all trimming steps (default: -1)"
)
parser.add_argument(
    "--orientation",
    metavar="int",
    default=0,
    type=int,
    help="re-orient reads (0: generic (default), 1: sense, -1: antisense)"
)

# output options
parser.add_argument(
    "--output_fusion",
    metavar="PATH",
    type=str,
    help="output fusion/chimeric reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--output_truncated",
    metavar="PATH",
    type=str,
    help="output truncated/non-full-length reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--output_full_length",
    metavar="PATH",
    type=str,
    help="output full-length reads to this file (use '-' for stdout)"
)
parser.add_argument(
    "--suffix_filtered",
    metavar="str",
    default=None,
    type=str,
    help="output filtered reads with the suffix"
)
