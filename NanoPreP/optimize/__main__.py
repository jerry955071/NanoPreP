from NanoPreP.optimize.argParser import parser
from NanoPreP.preptools.Annotator import Optimizer
from NanoPreP.seqtools.FastqIO import FastqIO, SeqFastq
from collections import Counter
from pathlib import Path
import plotly.graph_objects as go
import plotly.io
import gzip
import sys
import json

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
    params = parser.parse_args()
    
    # update params by user-provided config
    if params.config:
        config = json.load(open(params.config))
        params.update(config)
    
    # initiate Annotator object
    optimizer = Optimizer(
        p5_sense = params.p5_sense,
        p3_sense = params.p3_sense,
        isl5 = params.isl5,
        isl3 = params.isl3,
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
    with openg(params.input_fq, "r") as handle_in:
        # Sample at most `n` reads
        SAMPLED = 0
        for read in FastqIO.read(handle_in):
            # skip too-short reads if `skip_short`
            if params.skip_short > len(read):
                continue
            
            # skip low-quality reads if `skip_lowq` 
            if params.skip_lowq > SeqFastq.meanq(read):
                continue
            
            
            SAMPLED += 1
            if SAMPLED > params.n:
                SAMPLED -= 1
                break
            
            # primer alignment
            for _, (pid_isl, pid_rb) in optimizer.alignInB(read).items():
                counters[type1][pid_isl] += 1
                counters[type2][pid_rb] += 1
        
        # Plot sampling results
        pids = set(list(counters[type1].keys())) | set(list(counters[type2].keys()))
        pids = list(pids)
        pids.sort(reverse=True)
        for pid in pids:
            counters[type1][pid] += 0
            counters[type2][pid] += 0
            
        # Create figure
        fig = go.Figure()
        
        # Add trace `type2` data
        fig.add_trace(go.Bar(
            x=[f"{i:.2f}" for i in counters[type2].keys()],
            y=[v for k, v in counters[type2].items()],
            name=f"{type2} primer",
            marker_color=palatte[type2],
            opacity=.8,
            legendrank=2
        ))
        
        # Add trace `type1` data
        fig.add_trace(go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in counters[type1].keys()],
            y=[v for k, v in counters[type1].items()],
            name=f"{type1} primer",
            marker_color=palatte[type1],
            opacity=.8,
            legendrank=1
        ))

        # legendrank = 3
        hvtexts = []
        for pid in pids:
            tp = sum([v if k >= pid else 0 for k, v in counters[type1].items()])
            fp = sum([v if k >= pid else 0 for k, v in counters[type2].items()])
            precision = round(tp / (tp + fp), 3)
            hvtexts.append(f"True positives: {tp}<br>False positives: {fp}<br>Precision: {precision}")
            # fig.add_trace(
            #     go.Bar(
            #         x=[f"{pid:.2f}"],
            #         y=[max([v for k, v in counters[type1].items()] + [v for k, v in counters[type2].items()])],
            #         name=f"Cutoff: {round(pid, 2):.2f}",
            #         marker_color="red",
            #         hovertemplate=f"True positives: {tp}; False positives: {fp}; Precision: {precision}",
            #         width=.2,
            #         offset=-.5,
            #         visible="legendonly",
            #         legendrank=legendrank
            #     )
            # )
            # legendrank += 1
        fig.add_trace(go.Bar(
            visible=True,
            x=[f"{i:.2f}" for i in pids],
            y=[max([v for k, v in counters[type2].items()]+[v for k, v in counters[type1].items()])] * len(pids),
            marker_color="red",
            opacity=0,
            showlegend=False,
            hoverinfo="text",
            hovertext=hvtexts
        ))
        
        fig.update_xaxes(showspikes=True, spikecolor="red", spikesnap="cursor", spikemode="across")

        fig.update_layout(
            autosize=True,
            width=800,
            height=600,
            barmode='overlay', 
            xaxis_tickangle=-45,
            xaxis={'categoryorder':'category ascending'},
            font_family="Sans-serif",
            font_size=16,
            title={
                'text': f"Percent identity of {type1}/{type2} primer hits",
                'y':0.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            xaxis_title="Percent identity",
            yaxis_title="Counts",
            hovermode="x",
            hoverlabel=dict(namelength = -1)
        )
        
        with open(params.output, "w") as handle_html:
            out = """\
<html>
    <head><meta charset="utf-8" /></head>
    <style type="text/css">
        .gap-20 { 
            width:100%%; 
            height:20px; 
        }
        .gap-10 { 
            width:100%%; 
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
            font-size: 16px;
        }
    </style>
    <body>
        <div id="title">NanoPreP report</div>
        <div class="gap-20"></div>
        <div id="text">The figure below shows the percent identity of aligning primer sequences to ISL and read body on %(sampled) sampled reads generated with command: </div>
        <div class="gap-10"></div>
        <div id="code">nanoprep-optimizer %(args)</div>
        <div class="gap-10"></div>
        <div>%(plot)</div>
    </body>
</html>"""
            out = out.replace("%(args)", " ".join(sys.argv[1:]))
            out = out.replace("%(plot)", fig.to_html())
            out = out.replace("%(sampled)", f"{SAMPLED:,}")
            out = out.replace("%(type1)", type1)
            out = out.replace("%(type2)", type2)
            handle_html.write(out)
    return

if __name__ == "__main__":
    main()