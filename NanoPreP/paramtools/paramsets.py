Params = {
    "standard": {
        "skip_lowq": 7,
        "p5_sense": "TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG",
        "p3_sense": "A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC",
        "isl5": [0, 130],
        "isl3": [-60, -1],
        "pid_isl": 0.7,
        "pid_body": 0.7,
        "poly_w": 6,
        "poly_k": 4,
        "filter_short": 1,
        "trim_adapter": True,
        "trim_poly": True,
        "orientation": 1,
        "report": "report.json"
    },
    "annotate": {
        "p5_sense": "TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG",
        "p3_sense": "A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC",
        "isl5": [0, 130],
        "isl3": [-60, -1],
        "pid_isl": 0.7,
        "pid_body": 0.7,
        "poly_w": 6,
        "poly_k": 4,
        "orientation": 0,
        "output_full_length": "-",
        "output_fusion": "-",
        "output_truncated": "-",
        "report": "report.json"
    },
    "report": {
        "disable_annot": True,
        "report": "report.json"
    }
}
Defaults = {
    "skip_lowq": -1,
    "skip_short": -1,
    "disable_annot": False,
    "p5_sense": None,
    "p3_sense": None,
    "isl5": None,
    "isl3": None,
    "pid_isl": None,
    "pid_body": None,
    "poly_w": None,
    "poly_k": None,
    "orientation": 0,
    "trim_adapter": False,
    "trim_poly": False,
    "filter_short": -1,
    "filter_lowq": -1,
    "mode": None,
    "config": None,
    "report": None,
    "output_fusion": None,
    "output_truncated": None,
    "output_full_length": None,
    "suffix_filtered": None,
    "input_fq": None
}
