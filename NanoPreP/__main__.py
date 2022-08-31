#!usr/bin/env Python
from preptools.Annotator import Annotator
from preptools.Processor import Processor
from seqtools.FastqIO import FastqIO
from seqtools.SeqFastq import SeqFastq
from paramtools.paramsets import Params
from paramtools.argParser import parser
from datetime import datetime
import json


### parse arguments
args = parser.parse_args()

## collect parameters
params = {}

## get parameter presets
if args.mode:
    if args.mode in Params.keys():
        params = Params[args.mode]
    else:
        msg = "Available options to `mode`: "
        opts = ", ".join([i.__repr__() for i in Params.keys()])
        raise Exception(msg + opts)

## get parameter from config
if args.config:
    config = json.load(open(args.config))
    params.update(config)
    pass

## get command line arguments from ArgumentParser
for k, v in vars(args).items():
    if v:
        params[k] = v
    elif k not in params.keys():
        params[k] = v


## get time
## initiate report dict
report_dict = {
    "start time": datetime.now().strftime("%Y/%m/%d-%H:%M:%S"),
    "total reads": 0,
    "fusion/chimeric reads": 0,
    "truncated reads": 0,
    "full-length reads": 0,
    "too short": 0,
    "low quality": 0,
    "params": params,
    "stop time": None
}


# initiate Annotator
if not params["disable_annot"]:
    ## initiate Annotator
    annotator = Annotator(
        p5_sense=params["p5_sense"],
        p3_sense=params["p3_sense"],
        isl5=params["isl5"],
        isl3=params["isl3"],
        pid_isl=params["pid_isl"],
        pid_body=params["pid_body"],
        w=params["poly_w"],
        k=params["poly_k"]
    )

# open output files
handle_filtered = handle_passed = None
if params["output_filtered"]:
    open(params["output_filtered"], "w").close()
    handle_filtered = open(params["output_filtered"], "a") 

if params["output_passed"]:
    open(params["output_passed"], "w").close()
    handle_passed = open(params["output_passed"], "a")

# open `input_fq`
with open(params["input_fq"], "r") as handle_in:
    # stream processing
    for read in FastqIO.read(handle_in):
        PASS = True
        # annotate reads
        if not params["disable_annot"]:
            annotator.annotate(read)

        # fusion/chimeric reads
        if read.annot.fusion:
            if params["filter_fusion"]:
                PASS = False

        # full-length + truncated reads
        else:
            # try trimming
            if params["trim_poly"]:
                Processor.trimmer(read, True, True)
            elif params["trim_adapter"]:
                Processor.trimmer(read, False, True)

            # orient read
            if params["orientation"] != 0:
                Processor.orientor(read, to=params["orientation"])

            # check length
            if params["filter_short"] > len(read):
                PASS = False
                report_dict["too short"] += 1

            # check quality
            if params["filter_lowq"] > SeqFastq.meanq(read):
                PASS = False
                report_dict["low quality"] += 1

            # check full-length
            if params["filter_truncated"] and not read.annot.full_length:
                PASS = False


        # write to files
        if not PASS and handle_filtered:
            FastqIO.write(handle_filtered, read)

        if PASS and handle_passed:
            FastqIO.write(handle_passed, read)


        # update report_dict
        report_dict["total reads"] += 1
        report_dict["fusion/chimeric reads"] += read.annot.fusion
        report_dict["truncated reads"] += int(
            not read.annot.fusion and not read.annot.full_length
        )
        report_dict["full-length reads"] += read.annot.full_length



# get stopping time 
report_dict["stop time"] = datetime.now().strftime("%Y/%m/%d-%H:%M:%S")


# output report.json
with open(params["report"], "w") as handle:
    handle.write(json.dumps(report_dict, indent=4))
