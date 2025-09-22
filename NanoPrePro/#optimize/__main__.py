from NanoPrePro.optimize.argParser import parser, Defaults
from NanoPrePro.preptools.Annotator import Optimizer
from NanoPrePro.seqtools.FastqIO import FastqIO, SeqFastq
from NanoPrePro.optimize.template import template as html_template
from collections import Counter
from pathlib import Path
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import gzip
import sys
import json
import os

# create txt/gzip file handle
def openg(p:Path, mode:str):
    if isinstance(p, str):
        p = Path(p)
    if p.suffix == ".gz":
        return gzip.open(p, mode + "t")
    else:
        return open(p, mode)
    
    
# get primer identity distribution
def get_pid_counts(
        p5_sense:str,
        p3_sense:str,
        isl5:tuple,
        isl3:tuple,
        input_fq:str,
        n:int,
        skip_short:int,
        skip_lowq:int
    ):
    # initiate Annotator object
    optimizer = Optimizer(
        p5_sense =p5_sense,
        p3_sense = p3_sense,
        isl5 = isl5,
        isl3 = isl3,
    )
    
    # initiate counters
    counters = {
        "ISL": Counter(),
        "Read body": Counter()
    }
    
    # iterate through reads
    with openg(input_fq, "r") as handle_in:
        # Sample at most `n` reads
        SAMPLED = 0
        for read in FastqIO.read(handle_in):
            if SAMPLED > n:
                break
            
            # skip too-short reads if `skip_short`
            if skip_short > len(read):
                continue
            
            # skip low-quality reads if `skip_lowq` 
            if skip_lowq > SeqFastq.meanq(read):
                continue
                        
            # primer alignment
            for _, (pid_isl, pid_rb) in optimizer.alignInB(read).items():
                counters["ISL"][round(pid_isl, 2)] += 1
                counters["Read body"][round(pid_rb, 2)] += 1
                
            SAMPLED += 1
        

    # make sure counters share the same key list
    pids = set(list(counters["ISL"].keys())) | set(list(counters["Read body"].keys()))
    pids = list(pids)
    for pid in pids:
        counters["ISL"][pid] += 0
        counters["Read body"][pid] += 0
        
    return SAMPLED, counters


def plot_pid_counts(SAMPLED, counters, palatte, output, params, html_template):
    # Create figure
    fig = make_subplots(
        rows=1,
        cols=1,
        subplot_titles=[f"Sampled {SAMPLED:,} reads"]
    )
        
    # Add trace `Read body`
    fig.add_trace(
        go.Bar(
            x=[f"{i:.2f}" for i in counters["Read body"].keys()],
            y=[v for k, v in counters["Read body"].items()],
            name="Read body",
            marker_color=palatte["Read body"],
            opacity=.8,
            legendrank=2
        ),
        row=1, col=1
    )
        
    # Add trace `ISL`
    fig.add_trace(
        go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in counters["ISL"].keys()],
            y=[v for k, v in counters["ISL"].items()],
            name="ISL",
            marker_color=palatte["ISL"],
            opacity=.8,
            legendrank=1
        ),
        row=1, col=1
    )
        

    # Add trace `Precision` (opacity=0) and hovertext
    # make hovertexts
    pids = list(counters["ISL"].keys())
    hvtexts = []
    for pid in pids:
        tp = sum([v if k >= pid else 0 for k, v in counters["ISL"].items()])
        fp = sum([v if k >= pid else 0 for k, v in counters["Read body"].items()])
        precision = round(tp / (tp + fp), 3)
        hvtexts.append(
            f"True positives: {tp}<br>" +
            f"False positives: {fp}<br>" +
            f"Precision: {precision}"
        )
    # add trace
    fig.add_trace(
        go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in pids],
            y=[max([
                max(fig["data"][0]["y"]),
                max(fig["data"][1]["y"])
            ])] * len(pids),
            marker_color="red",
            opacity=0,
            showlegend=False,
            hoverinfo="text",
            hovertext=hvtexts
        ),
        row=1, col=1
    )
    
    
    # add spike
    fig.update_xaxes(
        showspikes=True,
        spikecolor="red",
        spikesnap="cursor",
        spikemode="across"
    )
    

    # update figure layout
    fig.update_layout(
        autosize=True,
        width=600,
        height=450,
        barmode='overlay', 
        xaxis={
            "categoryorder":"category ascending",
            "tickangle": -45,
            "title": {"text": "Percent identity"}
        },
        font_family="Sans-serif",
        font_size=16,
        yaxis_title="Counts",
        hovermode="x",
        hoverlabel=dict(namelength = -1)
    )
    
    
    # makedir for output
    os.makedirs(Path(output).parent, exist_ok=True)
    
    # write to file
    with open(output, "w") as handle_html:
        out = html_template
        out = out.replace("%(args)", " ".join(sys.argv[1:]))
        out = out.replace("%(plot)", fig.to_html())
        out = out.replace("%(params)", json.dumps(params))
        handle_html.write(out)
        
    return


def get_pid_cutoff(counters, precision_cutoff):
    # get pid cutoff
    pid_cutoff = None
    
    # iter over pids
    pids = list(counters["ISL"].keys())
    pids.sort(reverse=True)
    for pid in pids:
        tp = sum([v if k >= pid else 0 for k, v in counters["ISL"].items()])
        fp = sum([v if k >= pid else 0 for k, v in counters["Read body"].items()])
        precision = round(tp / (tp + fp), 3)
        if precision >= precision_cutoff:
            pid_cutoff = pid
        else:
            break

    return pid_cutoff


def main():    
    # load default parameters
    params = Defaults.copy()
    
    # parse arguments
    args = parser.parse_args()
    
    # update parameter from config
    if args.config:
        config = json.load(open(args.config))
        params.update(config)
        
    # update parameter with command line arguments
    for k, v in vars(args).items():
        if v != None:
            params[k] = v
            
    
    # get pid counts
    SAMPLED, counters = get_pid_counts(
        p5_sense=params["p5_sense"],
        p3_sense=params["p3_sense"],
        isl5=params["isl5"],
        isl3=params["isl3"],
        input_fq=params["input_fq"],
        n=params["n"],
        skip_short=params["skip_short"],
        skip_lowq=params["skip_lowq"]
    )


    # plot
    plot_pid_counts(
        SAMPLED=SAMPLED,
        counters=counters,
        palatte={
            "ISL": "#FFAF33",
            "Read body": "#2190FF"
        },
        output=params["output"],
        params=params,
        html_template=html_template
    )
    
if __name__ == "__main__":
    main()
