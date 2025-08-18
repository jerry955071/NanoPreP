# NanoPreP: a fully-equipped, fast, and memory-efficient pre-processor for ONT transcriptomic data

## Requirements
* Python (>= 3.7)  
* edlib (>=1.3.8)


## Getting started
**Option 1:** run with Docker
```
docker run chiaenu/nanoprep:latest nanoprep --help
```
**Option 2:** use pip install
```
pip install nanoprep-ffm
nanoprep --help
```

## General usage

NanoPreP optimizes adapter/primer identification parameters for each input file. During the optimization process, NanoPreP search for the best combination of (1) the adapter/primer substring used for alignment, (2) sequence similarity cutoff, and (3) aligned location cutoff that achieves the highest $F_{\beta}$ score in distinguishing real adapter/primer alignments and random alignments. 

The $\beta$ value in the formula of $F_{\beta}$ score greatly affect NanoPreP's behavior. The recommended range of $\beta$ is from 0.1 to 0.3, where smaller beta value lowers the chance of random alignments.


The general usage of NanoPreP to get ***high-quality, non-chimeric, full-length, strand-reoriented, adapter/primer-removed, polyA-removed*** reads:
```
nanoprep \
  --input_fq input.fq
  --beta 0.1 \
  --p5_sense 5_PRIMER_SEQUENCE \
  --p3_sense A{100}3_PRIMER_SEQUENCE \
  --trim_adapter \
  --trim_poly \
  --output_full_length output.fq \
  --report report.json
```
- `--input_fq input.fq` ← file contains raw sequences
- `--beta 0.1` ← optimize adapter/primer identification parameters using $F_{0.1}$ score
- `--p5_sense 5_PRIMER_SEQUENCE` ← 5' primer sequence in sense strand direction
- `--p3_sense A{100}3_PRIMER_SEQUENCE` ← expected length of polyA + 3' primer sequence in sense strand direction (see section [How to specify adapter/primer and polyA/T sequences](#HOWTO))
- `--output_full_length output.fq` ← write full-length reads to `output.fq`  
- `--report report.json` ← write details of the run to `report.json`  



<!-- TODO: why annotate reads? re-usable, time-saving, transparency, flexibility -->
After running this command, two output files `output.fq` and `report.json` will be written to your working directory.

The `report.json` records start/stop times, the parameters used, and the detail information of the input FASTQ file.  

The `output.fq` contains full-length reads processed by NanoPreP. For each processed read, NanoPreP appends the information of the read to the ID line (the line started with @): 
```
@read_1 strand=0.91 full_length=1 fusion=0 ploc5=0 ploc3=0 poly5=0 poly3=-20
AGAGGCTGGCGGGAACGGGC......TTTCAAAGCCAGGCGGATTC
+
+,),+'$)'%671*%('&$%......((&'(*($%$&%&$-((84*
```
As shown in the example above, several flags are used for the describe a read: 
|flag|regex|default|explanation|
|:-|-|-|:-|
|`strand`|-?\d+\.\d*|0|0: unknown; > 0: sense; < 0: antisense|
|`full_length`|[0\|1]|0|0: non-full-length; 1: full-length|
|`fusion`|[0\|1]|0|0: non-chimeric/-fusion; 1: chimeric/fusion|
|`ploc5`|-?\d+|-1|-1: unknown; 0: removed; > 0: 5' adapter/primer location|
|`ploc3`|-?\d+|-1|-1: unknown; 0: removed; > 0: 3' adapter/primer location|
|`poly5`|-?\d+|-1|0: unknown; > 0: 5' polymer length; < 0: trimmed 5' polymer length|
|`poly3`|-?\d+|-1|0: unknown; > 0: 3' polymer length; < 0: trimmed 3' polymer length|

According to the flags, the example "read1" is a sense strand (`strand=0.91`), full-length (`full_length=1`), non-chimeric (`fusion=0`),  adapter/primer removed (`ploc5=0 ploc3=0`), and polyA removed (`poly3=-20`) read.


## How to specify adapter/primer and polyA/T sequences <a id="HOWTO"></a>
Users need to provide the **adapter/primer** (and **polyA/T**) sequences to be searched for using options `--p5_sense` and `--p3_sense`. 

For example, the following command means that the 5' and 3' adatper/primer sequences on the sense strand are 'CATTC' and 'GACTA', respectively.
```
--p5_sense CATTC --p3_sense GACTA
```
If users wish to detect polyA/T tails, a pattern `N{M}` can be used to specify the location and length of polyA/T tails. The command below tells NanoPreP that there are poly`"A"` tails of a maximum length of `"50"` bases next to the 3' adapters/primers.
```
--p5_sense CATTC --p3_sense A{50}GACTA
```



## Full usage
```
usage: nanoprep [-h] [--version] --input_fq str [--config str] [--report str] [--processes int] [--batch_size int] [--seed int] [-n int] [--beta float] [--disable_annot] [--skip_lowq float] [--skip_short int] [--p5_sense str] [--p3_sense str]
                [--isl5 int int] [--isl3 int int] [--pid5 float] [--pid3 float] [--pid_body float] [--poly_w int] [--poly_k int] [--keep_adapter] [--keep_poly] [--filter_lowq float] [--filter_short int] [--orientation int] [--output_fusion str]
                [--output_truncated str] [--output_full_length str] [--suffix_filtered str]

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --input_fq str        input FASTQ
  --config str          use the parameters in this config file (JSON)
  --report str          output report file (JSON)
  --processes int       number of processes to use (default: 16)
  --batch_size int      number of records in each batch (default: 1000000)
  --seed int            seed for random number generator (default: 42)
  -n int                max number of reads to sample during optimization (default: 100000)
  --beta float          the beta parameter for the optimization (default: .1)
  --skip_lowq float     skip low-quality reads (default: 7)
  --skip_short int      skip too-short reads (default: 0)
  --p5_sense str        5' sense adatper/primer + polyA sequences
  --p3_sense str        3' sense adatper/primer + polyA sequences
  --isl5 int int        ideal searching location for 5' adapter/primer sequences (default: optimized)
  --isl3 int int        ideal searching location for 3' adapter/primer sequences (default: optimized)
  --pid5 float          5' adapter/primer percent identity cutoff (default: optimized)
  --pid3 float          3' adapter/primer percent identity cutoff (default: optimized)
  --pid_body float      adapter/primer percent identity cutoff (default: optimized)
  --poly_w int          window size for polyA/T identification
  --poly_k int          number of A/T to be expected in the window
  --trim_adapter        use this flag to trim adatper/primer sequences
  --trim_poly           use this flag to trim polyA/T sequences
  --filter_lowq float   filter low-quality reads after all trimming steps (default: 7)
  --filter_short int    filter too short reads after all trimming steps (default: 0)
  --orientation int     re-orient reads (0: generic , 1: sense (default), -1: antisense)
  --output_fusion str   output fusion/chimeric reads to this file (use '-' for stdout)
  --output_truncated str
                        output truncated/non-full-length reads to this file (use '-' for stdout)
  --output_full_length str
                        output full-length reads to this file (use '-' for stdout)
  --suffix_filtered str
                        output filtered reads with the suffix
```
