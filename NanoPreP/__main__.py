#!usr/bin/env Python
from pathlib import Path
from preptools.Annotator import Annotator
from preptools.Processor import Processor
from seqtools.FastqIO import FastqIO
from seqtools.SeqFastq import SeqFastq
from paramtools.paramsets import Params
from paramtools.argParser import parser
from datetime import datetime
import json


# parse arguments
args = parser.parse_args()

# collect parameters
params = {}

# get parameter presets
if args.mode:
    if args.mode in Params.keys():
        params = Params[args.mode]
    else:
        msg = "Available options to `mode`: "
        opts = ", ".join([i.__repr__() for i in Params.keys()])
        raise Exception(msg + opts)

# get parameter from config
if args.config:
    config = json.load(open(args.config))
    params.update(config)
    pass

# get command line arguments from ArgumentParser
for k, v in vars(args).items():
    if v:
        params[k] = v
    elif k not in params.keys():
        params[k] = v

# 
if not params["discard_short"]:
    params["discard_short"] = -1
if not params["discard_lowq"]:
    params["discard_lowq"] = -1
if not params["filter_short"]:
    params["filter_short"] = -1
if not params["filter_lowq"]:
    params["filter_lowq"] = -1



# initiate report dict
report_dict = {
    "start time": datetime.now().strftime("%Y/%m/%d-%H:%M:%S"),
    "total reads": 0,
    "fusion": {
        "passed": 0,
        "filtered": 0
    },
    "truncated": {
        "passed": 0,
        "filtered": 0
    },
    "full-length": {
        "passed": 0,
        "filtered": 0
    },
    "params": params,
    "stop time": None
}


# initiate an Annotator()
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
handle_out = {}
for name in ["output_fusion", "output_truncated", "output_full_length"]:
    for suff in ["suffix_passed", "suffix_filtered"]:
        cls = {
            "output_fusion": "fusion",
            "output_truncated": "truncated",
            "output_full_length": "full-length"
        }[name]
        passed = {
            "suffix_passed": "passed",
            "suffix_filtered": "filtered"
        }[suff]
        if params[name] and params[suff]:
            # open new file (clear if already exist)
            fout = Path(params[name])
            fout = fout.stem + "_" + params[suff] + fout.suffix
            open(fout, "w").close()
            # record to `handle_out`
            handle_out[(cls, passed)] = open(fout, "a")
        else:
            handle_out[(cls, passed)] = None


# open `input_fq`
with open(params["input_fq"], "r") as handle_in:
    # stream processing
    for read in FastqIO.read(handle_in):
        # add read count
        report_dict["total reads"] += 1

        # skip low-quality reads if `discard_lowq` 
        if params["discard_lowq"] > SeqFastq.meanq(read):
            continue

        PASS = "passed"
        # annotate reads
        if not params["disable_annot"]:
            annotator.annotate(read)

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
            PASS = "filtered"

        # check quality
        if params["filter_lowq"] > SeqFastq.meanq(read):
            PASS = "filtered"
        
        # classify reads into one of: fusion/full-length/truncated
        CLASS = None
        if read.annot.fusion:
            CLASS = "fusion"
        elif read.annot.full_length:
            CLASS = "full-length"
        else:
            CLASS = "truncated"
        
        # write to file if output file is assign
        if handle_out[(CLASS, PASS)]:
            FastqIO.write(handle_out[(CLASS, PASS)], read)

        # update report_dict
        report_dict[CLASS][PASS] += 1


# get the stopping time 
report_dict["stop time"] = datetime.now().strftime("%Y/%m/%d-%H:%M:%S")

# output report.json
with open(params["report"], "w") as handle:
    handle.write(json.dumps(report_dict, indent=4))
