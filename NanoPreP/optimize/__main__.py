from NanoPreP.optimize.argParser import parser, Defaults
from NanoPreP.preptools.Annotator import Optimizer
from NanoPreP.seqtools.FastqIO import FastqIO, SeqFastq
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

def main():
    # parse arguments
    args = parser.parse_args()
    
    # collect parameters
    params = Defaults.copy()
    
    # get parameter from config
    if args.config:
        config = json.load(open(args.config))
        params.update(config)
        
    # get command line arguments from ArgumentParser
    for k, v in vars(args).items():
        if v:
            params[k] = v

    
    # initiate Annotator object
    optimizer = Optimizer(
        p5_sense = params["p5_sense"],
        p3_sense = params["p3_sense"],
        isl5 = params["isl5"],
        isl3 = params["isl3"],
    )
    type1, type2 = "ISL", "Read body"
    counters = {
        type1: Counter(),
        type2: Counter()
    }
    
    
    # Color palette
    palatte = {
        type1: "#FFAF33",
        type2: "#2190FF"
    }
    
    
    #
    with openg(params["input_fq"], "r") as handle_in:
        # Sample at most `n` reads
        SAMPLED = 0
        for read in FastqIO.read(handle_in):
            # skip too-short reads if `skip_short`
            if params["skip_short"] > len(read):
                continue
            
            # skip low-quality reads if `skip_lowq` 
            if params["skip_lowq"] > SeqFastq.meanq(read):
                continue
            
            
            SAMPLED += 1
            if SAMPLED > params["n"]:
                SAMPLED -= 1
                break
            
            # primer alignment
            for _, (pid_isl, pid_rb) in optimizer.alignInB(read).items():
                counters[type1][round(pid_isl, 2)] += 1
                counters[type2][round(pid_rb, 2)] += 1
        

    # Plot sampling results
    pids = set(list(counters[type1].keys())) | set(list(counters[type2].keys()))
    pids = list(pids)
    for pid in pids:
        counters[type1][pid] += 0
        counters[type2][pid] += 0
        
    # Create figure
    fig = make_subplots(
        rows=1,
        cols=1,
        subplot_titles=[f"Sampled {SAMPLED:,} reads"]
    )
    
    # Add trace `type2` data to subplot (1, 1)
    fig.add_trace(
        go.Bar(
            x=[f"{i:.2f}" for i in counters[type2].keys()],
            y=[v for k, v in counters[type2].items()],
            name=f"{type2}",
            marker_color=palatte[type2],
            opacity=.8,
            legendrank=2
        ),
        row=1, col=1
    )
    
    # Add trace `type1` data
    fig.add_trace(
        go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in counters[type1].keys()],
            y=[v for k, v in counters[type1].items()],
            name=f"{type1}",
            marker_color=palatte[type1],
            opacity=.8,
            legendrank=1
        ),
        row=1, col=1
    )
    

    # legendrank = 3
    hvtexts = []
    for pid in pids:
        tp = sum([v if k >= pid else 0 for k, v in counters[type1].items()])
        fp = sum([v if k >= pid else 0 for k, v in counters[type2].items()])
        precision = round(tp / (tp + fp), 3)
        hvtexts.append(
            f"True positives: {tp}<br>" +
            f"False positives: {fp}<br>" +
            f"Precision: {precision}"
        )
        
    fig.add_trace(
        go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in pids],
            y=[max([
                max(fig["data"][0]["y"]),
                max(fig["data"][0]["y"])
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
    
    
    # makedir
    os.makedirs(Path(params["output"]).parent, exist_ok=True)
    with open(params["output"], "w") as handle_html:
        out = """\
<html>
    <head><meta charset="utf-8" /></head>
    <style type="text/css">
        .gap-20 { 
            width:100%; 
            height:20px; 
        }
        .gap-10 { 
            width:100%; 
            height:10px; 
        } 
        #title {
            font-family: Sans-serif;
            font-size: 24px;
            font-weight: 600;
        }
        #text {
            font-family: Sans-serif;
            font-size: 16px;
            font-weight: 600;
        }
        #comment {
            font-family: Sans-serif;
            font-size: 12px;
        }
        #code {
            font-family: Monospace;
            font-size: 14px;
        }
    </style>
    <body>
        <div id="title">NanoPreP-optimize</div>
        <div class="gap-20"></div>
        <div>%(plot)</div>
        <div id="text">Command line:</div>
        <div id="code">nanoprep-optimize %(args)</div>
        <div class="gap-10"></div>
        <div id="text">Parameters:</div>
        <div id="code">%(params)</div>
    </body>
</html>"""
        out = out.replace("%(args)", " ".join(sys.argv[1:]))
        out = out.replace("%(plot)", fig.to_html())
        out = out.replace("%(params)", json.dumps(params))
        handle_html.write(out)
    return

if __name__ == "__main__":
    main()