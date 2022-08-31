from argparse import ArgumentParser

# initiate ArgumentParser
parser = ArgumentParser(description='Arguments availible to NanoPreP')

# annotaion options
parser.add_argument(
    "--disable_annot",
    action="store_true",
    help="disable annotation for NanoPreP-annotated FASTQ)"
)
parser.add_argument(
    "--p5_sense",
    metavar="",
    type=str,
    help="5' sense adatper/primer/polymer sequences"
)
parser.add_argument(
    "--p3_sense",
    metavar="",
    type=str,
    help="3' sense adatper/primer/polymer sequences"
)
parser.add_argument(
    "--isl5",
    metavar="",
    nargs=2,
    type=int,
    help="ideal searching location for 5' adapter/primer sequences"
)
parser.add_argument(
    "--isl3",
    metavar="",
    nargs=2,
    type=int,
    help="ideal searching location for 3' adapter/primer sequences"
)
parser.add_argument(
    "--pid_isl",
    metavar="",
    type=float,
    help="adapter/primer percent identity cutoff (in ISL)"
)
parser.add_argument(
    "--pid_body",
    metavar="",
    type=float,
    help="adapter/primer percent identity cutoff (on read body)"
)
parser.add_argument(
    "--poly_w",
    metavar="",
    type=int,
    help="window size for homopolymer identification"
)
parser.add_argument(
    "--poly_k",
    metavar="",
    type=int,
    help="number of monomers to be expected in the window"
)

# processing (filtering/trimming/orientation) options
parser.add_argument(
    "--filter_short",
    metavar="",
    default=0,
    type=int,
    help="filter too short reads after all trimming steps (default: 0)"
)
parser.add_argument(
    "--filter_lowq",
    metavar="",
    default=-1,
    type=int,
    help="filter low-quality reads after all trimming steps (default: -1)"
)
parser.add_argument(
    "--filter_fusion",
    action="store_true",
    help="filter fusion/chimeric reads"
)
parser.add_argument(
    "--filter_truncated",
    action="store_true",
    help="filter truncated/non-full-length reads"
)
parser.add_argument(
    "--trim_adapter",
    action="store_true",
    help="trim adatper/primer sequences"
)
parser.add_argument(
    "--trim_poly",
    action="store_true",
    help="trim homopolymers"
)
parser.add_argument(
    "--orientation",
    metavar="",
    default=0,
    help="re-orient reads (0: generic (default), 1: sense, -1: antisense)"
)

# general options
parser.add_argument(
    "input_fq",
    metavar="input.fq",
    type=str,
    help="input FASTQ (required)"
)
parser.add_argument(
    "--output_filtered",
    metavar="",
    type=str,
    help="output filter reads FASTQ"
)
parser.add_argument(
    "--output_passed",
    metavar="",
    type=str,
    help="output passed reads FASTQ"
)
parser.add_argument(
    "--report",
    metavar="",
    type=str,
    help="output report JSON file"
)
parser.add_argument(
    "--config",
    metavar="",
    type=str,
    help="provide parameters in JSON format "
    "(can be overriden by command line arguments)"
)
parser.add_argument(
    "--mode",
    metavar="",
    type=str,
    help="presets of processing parameters "
    "[strandard|untrimmed|polyA|annotate] "
    "(can be overriden by `config` and command line arguments) "
)
