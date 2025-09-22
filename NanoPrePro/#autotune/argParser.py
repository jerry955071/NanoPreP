from argparse import ArgumentParser
# initiate ArgumentParser
parser = ArgumentParser(description='Arguments availible to nanoprep-autotune')

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
    default=10000,
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
    "--output_report",
    metavar="PATH",
    default="report.html",
    type=str,
    help="output report.html (default: report.html)"
)
parser.add_argument(
    "--output_config",
    metavar="PATH",
    default="config.json",
    type=str,
    help="output config.json (default: config.json)"
)

# annotation options
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
    required=True,
    help="5' sense adatper/primer sequences"
)
parser.add_argument(
    "--p3_sense",
    metavar="str",
    type=str,
    required=True,
    help="3' sense adatper/primer sequences"
)

