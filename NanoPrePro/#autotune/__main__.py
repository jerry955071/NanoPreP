from NanoPrePro.autotune.argParser import parser
from NanoPrePro.aligntools.edlibAligner import edlibAligner as aligner
from NanoPrePro.preptools.Annotator import Tuner
from NanoPrePro.seqtools.FastqIO import FastqIO
from NanoPrePro.seqtools.SeqFastq import SeqFastq
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import gzip

# open files
def openg(p:Path, mode:str):
    if isinstance(p, str):
        p = Path(p)
    if p.suffix == ".gz":
        return gzip.open(p, mode + "t")
    else:
        return open(p, mode)

# main
def main():
    # parse arguments
    params = parser.parse_args()
    
    # update params by user-provided config
    if params.config:
        config = json.load(open(params.config))
        params.update(config)
    
    # update params of antisense primers
    params.p5_anti = SeqFastq.reverse_complement_static(params.p3_sense)
    params.p3_anti = SeqFastq.reverse_complement_static(params.p5_sense)
    
    # initiate annotator
    tuner = Tuner(
        params.p5_sense,
        params.p5_anti,
        params.p3_sense,
        params.p3_anti 
    )
    max_plen = max(
        len(params.p5_sense),
        len(params.p3_sense)
    )
    plens = [i for i in range(20, max_plen, 10)]
    plens[-1] = max_plen

    
    # assess primer to primer similarity
    p2pid = {plen: tuner.primer2primer(plen) for plen in plens}
    json.dump(p2pid, open("test/p2p_ident.json", "w"))    
    
    # find pid cutoff
    p2rid = {
        plen: {
            "p5": {"pid": [], "location": [], "type": []},
            "p3": {"pid": [], "location": [], "type": []} 
        } 
        for plen in plens
    }

    p2rid = {
        "rid": [],
        "end": [],
        "plen": [],
        "pid": [],
        "location": [],
        "type": []
    }
    SKIPPED = 0
    
    with openg(params.input_fq, "r") as handle_in:
        N = 0
        for read in FastqIO.read(handle_in):
            # count +1
            N += 1
            # check n
            if N > params.n:
                break
            
            # skip too-short reads if `skip_short`
            if params.skip_short > len(read):
                SKIPPED += 1
                continue

            # skip low-quality reads if `skip_lowq` 
            if params.skip_lowq > SeqFastq.meanq(read):
                SKIPPED += 1
                continue
            
            # tuning
            for plen in plens:
                res = tuner.primer2read(read, plen)
                p2rid["rid"] += [read.id] * 8
                p2rid["end"] += ["p5"] * 4 + ["p3"] * 4
                p2rid["plen"] += [plen] * 8
                p2rid["type"] += ["primary", "secondary", "primary-competing", "secondary-competing"] * 2
                # 5' end
                if res["p5_sense"][0].pid > res["p5_anti"][0].pid:
                    p5 = res["p5_sense"]
                    p5_c = res["p5_anti"]
                else:
                    p5 = res["p5_anti"]
                    p5_c = res["p5_sense"]
                # 3' end
                if res["p3_sense"][0].pid > res["p3_anti"][0].pid:
                    p3 = res["p3_sense"]
                    p3_c = res["p3_anti"]
                else:
                    p3 = res["p3_anti"]
                    p3_c = res["p3_sense"]
                    
                p2rid["pid"] += [p5[0].pid, p5[1].pid, p5_c[0].pid, p5_c[1].pid]
                p2rid["location"] += [p5[0].location, p5[1].location, p5_c[0].location, p5_c[1].location]
                p2rid["pid"] += [p3[0].pid, p3[1].pid, p3_c[0].pid, p3_c[1].pid]
                p2rid["location"] += [p3[0].location, p3[1].location, p3_c[0].location, p3_c[1].location]
                
                
        # plotting
        p2rid = pd.DataFrame(p2rid)
        p2rid["pid"] = p2rid["pid"].apply(np.round, False, args=[2])
        p2rid.to_csv("test/data.tsv", sep="\t", header=True, index=False)
        # site = "p5"
        # for plen in plens:
        #     fig, ax = plt.subplots(figsize=(8, 6))
        #     data = p2rid[plen][site]
        #     # sns.scatterplot(x="pid", y="location",
        #     #         hue="type",
        #     #         palette="ch:r=-.2,d=.3_r",
        #     #         s=1, linewidth=0,
        #     #         alpha=.9,
        #     #         data=data, ax=ax)
        #     plt.yscale("log")
        #     # sns.countplot(
        #     #     data=data,
        #     #     x="pid",
        #     #     hue="type",
        #     #     palette="ch:r=-.2,d=.3_r",
        #     #     ax=ax
        #     # )
        #     "linear', 'log', 'symlog', 'asinh', 'logit', 'function', 'functionlog'"
        #     # sns.histplot(
        #     #     data=data,
        #     #     x="location",
        #     #     hue="type",
        #     #     palette="ch:r=-.2,d=.3_r",
        #     #     ax=ax
        #     # )
        #     sns.stripplot(
        #         data=data,
        #         x="pid",
        #         y="location",
        #         hue="type",
        #         hue_order=["primary", "secondary"],
        #         palette="ch:r=-.2,d=.3_r",
        #         s=3
        #     )
        #     ax.set_xticklabels(
        #         [round(float(i.get_text()), 2) for i in ax.get_xticklabels()],
        #         rotation=70
        #     )
            
        #     fig.savefig(f"test/fig_{site}_loc2pid_{plen}.jpg")

        # print(0)
    
if __name__ == "__main__":
    main()