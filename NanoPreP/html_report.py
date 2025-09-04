from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import json, sys


class HTML_report:
    def __init__(self):
        self.template = """\
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
</html>
"""
    def write(self, path):
        with open(path, "w") as handle_html:
            handle_html.write(self.template)
        return
    
    def update_optimized_data(self, data, n_iqr, palette={"positive": "blue", "negative": "orange"}):
        # 1. Preprocess optimization data
        # one alignment one row --> record counts of identical alignment results
        data_agg = data.groupby(["pid","loc","cls","plen","site"]).size().reset_index(name="count")
        # add jitter
        data_agg["pid_jitter"] = data_agg["pid"] + np.random.uniform(-0.005,0.005,len(data_agg))

        # 2. Create plot (1 row, 2 col)
        fig = make_subplots(rows=1, cols=2, subplot_titles=["5' AP", "3' AP"], shared_yaxes=True)
        
        # 3. Add traces
        plens = sorted(data_agg["plen"].unique())
        sites = ["left", "right"]
        for plen_idx, plen in enumerate(plens):
            for site_idx, site in enumerate(sites):
                subset = data_agg[(data_agg["plen"]==plen) & (data_agg["site"]==site)]
                # plot positive/negative alingment results
                for cls in ["positive", "negative"]:
                    cls_subset = subset[subset["cls"]==cls]
                    fig.add_trace(
                        go.Scatter(
                            x=cls_subset["pid_jitter"],
                            y=cls_subset["loc"],
                            mode="markers",
                            marker=dict(
                                color=palette[cls],
                                size=np.clip(np.log1p(cls_subset["count"])*3,2,12),
                                line_width=0,
                                opacity=0.6
                            ),
                            name={"positive": "Real AP (1st alignment)", "negative": "Random AP (2nd alignment)"}[cls],  # legend entry
                            hoverinfo="skip",  # suppress hover for raw dots
                            visible=(plen_idx==0),
                            showlegend=(site == "left")
                        ),
                        row=1, col=site_idx+1
                    )
                # calculate metrics
                pids = subset["pid"].unique()
                tmp_loc = data[(data["cls"] == "positive") & (data["plen"] == plen) & (data["site"] == site)]["loc"]
                locs = np.linspace(tmp_loc.median(), tmp_loc.quantile(.9), n_iqr)
                metrics_rows = []
                for i in pids:
                    for j in locs:
                        accepted = (subset["pid"] >= i) & (subset["loc"] <= j)
                        is_positive = subset["cls"] == "positive"
                        tp = subset[is_positive & accepted]["count"].sum()
                        fp = subset[~is_positive & accepted]["count"].sum()
                        fn = subset[is_positive & ~accepted]["count"].sum()
                        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
                        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                        fscore = 2 * prec * recall / (prec + recall) if (prec + recall > 0) else 0
                        metrics_rows.append([i, j, prec, recall, fscore])
                metrics_df = pd.DataFrame(metrics_rows)
                
                # plot metrics data (invisible; for hovertext)
                fig.add_trace(
                    go.Scatter(
                        x=metrics_df[0] if not metrics_df.empty else [],
                        y=metrics_df[1] if not metrics_df.empty else [],
                        mode="markers",
                        marker=dict(opacity=0),
                        customdata=metrics_df.values if not metrics_df.empty else None,
                        hovertemplate=(
                            "Similarity cutoff: %{x}<br>Location cutoff: %{y}<br>"
                            "Precision: %{customdata[2]:.2f}<br>"
                            "Recall: %{customdata[3]:.2f}<br>"
                            "F1: %{customdata[4]:.2f}<br><extra></extra>"
                        ),
                        hoverlabel=dict(
                            align="left",
                            bgcolor="whitesmoke"
                        ),
                        visible=(plen_idx==0),
                        showlegend=False
                    ),
                    row=1, col=site_idx+1
                )
        
        # 4. Adjust layouts
        fig.update_layout(
            hovermode="closest",
            xaxis_title="Sequence similarity<br>(Percent identity)",
            xaxis2_title="Sequence similarity<br>(Percent identity)",
            yaxis_title="Location<br>(<i>n</i>-bp from read termini)",
            width=1000, height=500,
        )
        fig.update_yaxes(
            autorange="reversed",
            type="log",
            showspikes=True,
            spikecolor="red",
            spikemode="across",
            spikesnap="cursor"
        )
        fig.update_xaxes(
            showspikes=True,
            spikecolor="red",
            spikemode="across",
            spikesnap="cursor"
        )
        
        # 5. Add slider
        num_sites = len(sites)
        traces_per_site = 3  # positive, negative, metrics
        traces_per_plen = num_sites * traces_per_site
        total_traces = len(fig.data)
        steps = []
        for plen_idx, plen in enumerate(plens):
            visible = [False] * total_traces
            base = plen_idx * traces_per_plen
            for s in range(num_sites):
                for offset in range(traces_per_site):
                    idx = base + s*traces_per_site + offset
                    if idx < total_traces:
                        visible[idx] = True
            steps.append(dict(method="update", args=[{"visible": visible}], label=str(plen)))
        fig.update_layout(
            sliders=[dict(
                active=0,
                currentvalue={"prefix":"AP length (%): "},
                pad={"t":80},
                steps=steps,
                font=dict(size=14)
            )]
        )
        
        # 6. Update template
        self.template.replace("%(plot)", fig.to_html())
        
        return

    def update_command_line(self, cmdline):
        self.template.replace("%(args)", " ".join(sys.argv[1:]))
        return
    
    def update_params(self, params):
        self.template.replace("%(params)", json.dumps(params))
        return