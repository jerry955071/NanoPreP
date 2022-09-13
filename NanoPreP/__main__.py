#!usr/bin/env Python
from pathlib import Path
import re
from NanoPreP.preptools.Annotator import Annotator
from NanoPreP.preptools.Processor import Processor
from NanoPreP.seqtools.FastqIO import FastqIO
from NanoPreP.seqtools.SeqFastq import SeqFastq
from NanoPreP.paramtools.paramsets import Params, Defaults
from NanoPreP.paramtools.argParser import parser
from datetime import datetime
import sys
import json

def main():
    # parse arguments
    args = parser.parse_args()

    # collect parameters
    params = Defaults.copy()

    # get parameter presets
    if args.mode:
        if args.mode in Params.keys():
            params.update(Params[args.mode])
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
        if not v == Defaults[k]:
            params[k] = v    
        elif f"--{k}" in sys.argv:
            params[k] = v
        

    # initiate report dict
    report_dict = {
        "start time": datetime.now().strftime("%Y/%m/%d-%H:%M:%S"),
        "total reads": 0,
        "skipped": 0,
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
        "stop time": None,
        "params": params
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
        cls = {
            "output_fusion": "fusion",
            "output_truncated": "truncated",
            "output_full_length": "full-length"
        }[name]
        # output passed
        if params[name] == "-":
            handle_out[(cls, "passed")] = sys.stdout
        elif params[name]:
            # open new file (clear if already exist)
            fout = Path(params[name])
            open(fout, "w").close()
            handle_out[(cls, "passed")] = open(fout, "a")
        else:
            handle_out[(cls, "passed")] = None

        # output filtered
        if params[name] and params["suffix_filtered"]:
            # open new file (clear if already exist)
                fout = Path(params[name])
                fout = fout.stem + "_" + params["suffix_filtered"] + fout.suffix
                open(fout, "w").close()
                # record to `handle_out`
                handle_out[(cls, "filtered")] = open(fout, "a")
        else:
            handle_out[(cls, "filtered")] = None


    # open `input_fq`
    with open(params["input_fq"], "r") as handle_in:
        # stream processing
        for read in FastqIO.read(handle_in):
            # add read count
            report_dict["total reads"] += 1

            # skip too-short reads if `skip_short`
            if params["skip_short"] > len(read):
                report_dict["skipped"] += 1
                continue

            # skip low-quality reads if `skip_lowq` 
            if params["skip_lowq"] > SeqFastq.meanq(read):
                report_dict["skipped"] += 1
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
            
            # write to file if output specified
            if handle_out[(CLASS, PASS)]:
                FastqIO.write(handle_out[(CLASS, PASS)], read)

            # update report_dict
            report_dict[CLASS][PASS] += 1


    # get the stopping time 
    report_dict["stop time"] = datetime.now().strftime("%Y/%m/%d-%H:%M:%S")


    # output report.json
    if params["report"]:
        with open(params["report"], "w") as handle:
            handle.write(json.dumps(report_dict, indent=4))

if __name__ == "__main__":
    main()