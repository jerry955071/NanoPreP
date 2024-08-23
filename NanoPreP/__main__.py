from NanoPreP.preptools.Annotator import Annotator
from NanoPreP.preptools.Processor import Processor
from NanoPreP.seqtools.FastqIO import FastqIO, FastqIndexIO
from NanoPreP.seqtools.SeqFastq import SeqFastq
# from NanoPreP.paramtools.paramsets import Params, Defaults
from NanoPreP.paramtools.argParser import parser
from NanoPreP.preptools.Optimizer import Optimizer
from datetime import datetime
from pathlib import Path
import multiprocessing as mp
import os, sys, json, gzip, random, logging

# TODO: 1. test for the number of reads to sample for optimization

# initiate global variable `PARAMS`
PARAMS = {}

# set logging style
logging.basicConfig(
    format="[%(asctime)s] PID=%(process)d: %(message)s",
    datefmt="%H:%M:%S",
    level=logging.DEBUG
)


# main function
def main():
    # get parameters
    num_reads = get_params()
    
    # create output queue and process pool
    manager = mp.Manager()
    output_queue = manager.Queue()
    pool = mp.Pool(processes=PARAMS["processes"] + 1) # +1 for `fq_writer`
    
    # initiate `REPORT_DICT`
    REPORT_DICT = {}
    REPORT_DICT["start time"] = datetime.now().strftime("%Y/%m/%d-%H:%M:%S")
    REPORT_DICT["total reads"] = 0
    REPORT_DICT["skipped"] = 0
    REPORT_DICT["fusion/passed"] = 0
    REPORT_DICT["fusion/filtered"] = 0
    REPORT_DICT["truncated/passed"] = 0
    REPORT_DICT["truncated/filtered"] = 0
    REPORT_DICT["full-length/passed"] = 0
    REPORT_DICT["full-length/filtered"] = 0
    REPORT_DICT["stop time"] = None
    
    # put fq_writer/report_logger to work
    watcher_out = pool.apply_async(fq_writer, (output_queue,))

    # create batch tasks
    tasks = []
    batches = max(num_reads // PARAMS["batch_size"] + 1, PARAMS["processes"])
    batch_size = num_reads // batches + 1
    for batch_id in range(batches):
        tasks.append((
            PARAMS["input_fq"],
            output_queue,
            batch_id,
            batch_size
        ))

    # collect results from the workers through the pool result queue
    read_counts = pool.starmap_async(batch_worker, tasks)
    for read_count in read_counts.get():
        for k, v in read_count.items():
            REPORT_DICT[k] += v

    # close the pool and wait for the workers to finish
    output_queue.put('kill')
    pool.close()
    pool.join()
    
    # get stop time
    REPORT_DICT["stop time"] = datetime.now().strftime("%Y/%m/%d-%H:%M:%S")

    # write report.json
    if PARAMS["report"]:
        with openg(Path(PARAMS["report"]), "w") as handle:
            REPORT_DICT["params"] = PARAMS
            handle.write(json.dumps(REPORT_DICT, indent=4))
            
    return


# worker function
def batch_worker(input_file, output_queue, batch_id, batch_size):
    # logging
    logging.info("Start batch %d" % batch_id)
    
    # initiate `read_counter`
    read_counter = {
        "total reads": 0,
        "skipped": 0,
        "fusion/passed": 0,
        "fusion/filtered": 0,
        "truncated/passed": 0,
        "truncated/filtered": 0,
        "full-length/passed": 0,
        "full-length/filtered": 0
    }
    
    # initiate Annotator
    annotator = Annotator(
        p5_sense=PARAMS["p5_sense"],
        p3_sense=PARAMS["p3_sense"],
        isl5=PARAMS["isl5"],
        isl3=PARAMS["isl3"],
        pid5=PARAMS["pid5"],
        pid3=PARAMS["pid3"],
        pid_body=PARAMS["pid_body"],
        w=PARAMS["poly_w"],
        k=PARAMS["poly_k"]
    )
    
    # get reads (stored in memory)
    logging.info(f"Batch {batch_id}: Loading reads")
    reads = FastqIO.batch_read(
        input_file,
        batch_id * batch_size,
        batch_size
    )
    logging.info(f"Batch {batch_id}: Loaded {len(reads):,d} reads")
    
    # iterate through reads
    for read in reads:
        # update read counter
        read_counter["total reads"] += 1
        
        # read = indexed_fq.get(rname)
        # skip too-short reads if `skip_short`
        if PARAMS["skip_short"] > len(read):
            read_counter["skipped"] += 1
            continue

        # skip low-quality reads if `skip_lowq` 
        if PARAMS["skip_lowq"] > SeqFastq.meanq(read):
            read_counter["skipped"] += 1
            continue

        PASS = "passed"
        # annotate reads
        if not PARAMS["disable_annot"]:
            annotator.annotate(read)

        # try trimming
        if not PARAMS["keep_poly"]:
            Processor.trimmer(read, True, True)
        elif not PARAMS["keep_adapter"]:
            Processor.trimmer(read, False, True)

        # orient read
        if PARAMS["orientation"] != 0:
            Processor.orientor(read, to=PARAMS["orientation"])

        # check length
        if PARAMS["filter_short"] > len(read):
            PASS = "filtered"

        # check quality
        if PARAMS["filter_lowq"] > SeqFastq.meanq(read):
            PASS = "filtered"
        
        # classify reads into one of: fusion/full-length/truncated
        CLASS = None
        if read.annot.fusion:
            CLASS = "fusion"
        elif read.annot.full_length:
            CLASS = "full-length"
        else:
            CLASS = "truncated"
        
        # put read into queue
        output_queue.put([read, (CLASS, PASS)])
        
        # update read counter
        read_counter[f"{CLASS}/{PASS}"] += 1
    
    logging.info(f"Finished batch {batch_id}")
    return read_counter


# listener function
def fq_writer(queue):
    # initiate handle_dict
    handle_dict = {}
    
    # create/open output files
    for name in ["output_fusion", "output_truncated", "output_full_length"]:
        cls = {
            "output_fusion": "fusion",
            "output_truncated": "truncated",
            "output_full_length": "full-length"
        }[name]
        # output passed
        if PARAMS[name] == "-":
            handle_dict[(cls, "passed")] = sys.stdout
        elif PARAMS[name]:
            os.makedirs(Path(PARAMS[name]).parent, exist_ok=True)
            handle_dict[(cls, "passed")] = openg(Path(PARAMS[name]), "w")
        else:
            handle_dict[(cls, "passed")] = None

        # output filtered
        if PARAMS[name] and PARAMS["suffix_filtered"]:
            # open new file (truncate if already exist)
            fout = Path(PARAMS[name])
            fout = fout.stem + "_" + PARAMS["suffix_filtered"] + fout.suffix
            handle_dict[(cls, "filtered")] = openg(fout, "w")
        else:
            handle_dict[(cls, "filtered")] = None
    
    # write reads to output files
    while True:
        m = queue.get()
        if m == "kill":
            break
        read, (cls, pass_) = m
        if handle_dict[(cls, pass_)]:
            FastqIO.write(handle_dict[(cls, pass_)], read)
            
    # close handles
    for handle in handle_dict.values():
        handle.close()
        
    return


# get parameters from command line arguments
def get_params():
    global PARAMS
    
    # parse arguments
    args = parser.parse_args()
    PARAMS = vars(args)

    
    # update `params` from config (mostly for rerun)
    if args.config:
        config = json.load(open(args.config))
        PARAMS.update(config)


    # optimize AP identification parameters if `--beta` is specified   
    if PARAMS["beta"]:
        # initialize optimizer
        optimizer = Optimizer(
            p5_sense=PARAMS["p5_sense"],
            p3_sense=PARAMS["p3_sense"]
        )
        
        # sample `n` reads
        logging.info("Counting records in FASTQ file")
        sampled_fq, num_reads = FastqIO.sample(
            PARAMS["input_fq"],
            PARAMS["n"], 
            PARAMS["seed"]
        )
        logging.info(f"Found {num_reads:,d} records in FASTQ file")
        logging.info(f"Sampled {len(sampled_fq):,d} records for optimization")
    
        # optimize parameters            
        out = optimizer.optimize(
            fq_iter=sampled_fq,
            plens=[.5, .6, .7, .8, .9, 1.0],
            n_iqr=[.5, 1, 1.5, 2],
            processes=PARAMS["processes"],
            target="fscore",
            beta=PARAMS["beta"]
        )
        
        # update `PARAMS`
        PARAMS["p5_sense"] = out["left"]["seq"]
        PARAMS["p3_sense"] = out["right"]["seq"]
        PARAMS["pid5"] = out["left"]["pid"]
        PARAMS["pid3"] = out["right"]["pid"]
        PARAMS["pid_body"] = max(
            out["left"]["pid"],
            out["right"]["pid"]
        )
        PARAMS["isl5"] = (0, int(out["left"]["loc"]))
        PARAMS["isl3"] = (-int(out["right"]["loc"]), -1)
        
        res = (
            f"Optimization results:\n"
            f"5' primer:\n"
            f"  - sequence={PARAMS['p5_sense']}\n"
            f"  - percent identity cutoff={PARAMS['pid5']}\n"
            f"  - ideal searching location=({PARAMS['isl5'][0]}, {PARAMS['isl5'][1]})\n"
            f"  - fscore={out['right']['fscore']:.2f}\n"
            f"  - precision={out['left']['prec']:.2f}\n"
            f"  - recall={out['left']['recall']:.2f}\n"
            f"3' primer:\n"
            f"  - sequence={PARAMS['p3_sense']}\n"
            f"  - percent identity cutoff={PARAMS['pid3']}\n"
            f"  - ideal searching location=({PARAMS['isl3'][0]}, {PARAMS['isl3'][1]})\n"
            f"  - fscore={out['right']['fscore']:.3f}\n"
            f"  - precision={out['right']['prec']:.3f}\n"
            f"  - recall={out['right']['recall']:.3f}"
        )
        logging.info(res)
            
    return num_reads


# extended open function
def openg(p:Path, mode:str):
    if p.suffix == ".gz":
        return gzip.open(p, mode + "t")
    else:
        return open(p, mode)
    
if __name__ == "__main__":
   main()
   logging.info("Finished all batches.")