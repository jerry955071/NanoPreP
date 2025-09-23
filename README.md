# NanoPrePro: a fully-equipped, fast, and memory-efficient pre-processor for ONT transcriptomic data

## Requirements
* Python (>= 3.7)  
* edlib (>=1.3.8)
* Pandas
* NumPy
* Plotly

## Installation
**Option 1:** run with Docker
```
docker run chiaenu/nanoprepro:latest nanoprepro --help
```
**Option 2:** use pip install
```
pip install nanoprepro
nanoprepro --help
```

## Quick start

Run a standard preprocessing pipeline using NanoPrePro as follows:

```
nanoprepro \
  --input_fq input.fq \
  --beta 0.2 \
  --p5_sense 5_PRIMER_SEQUENCE \
  --p3_sense A{100}3_PRIMER_SEQUENCE \
  --filter_lowq 7 \
  --trim_adapter \
  --trim_poly \
  --orientation 1 \
  --output_full_length output.fq \
  --report report.html
```

This command performs the following preprocessing steps and generates a report file (`report.html`):

1. `--beta 0.2`: performs $F_{\beta=0.2}$ optimization for adapter/primer alignment cutoffs.
2. `--output_full_length output.fq`: identifies full-length reads and write to file.
3. `--trim_adapter`: trims adapter/primer sequences.
4. `--trim_poly`: trims poly(A/T) sequences.
5. `--orientation 1`: reorients reads.
6. `--filter_lowq 7`: filters low-quality reads.


## Documentation
For detailed usage, please refer to the [documentation](https://nanoprepro.readthedocs.io/en/latest/index.html).
