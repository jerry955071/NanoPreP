# NanoPreP
NanoPreP is a fully-functional, fast, memory-efficient pre-processor for ONT transcriptomic data
<br>  

## Requirements
* Python (>= 3.7)  
* edlib (1.3.8) for Python
<br>  

## Getting started
```
git clone https://github.com/Woodformation1136/NanoPreP.git
cd NanoPreP
python NanoPreP --help
```
<br>  

## Workflow
1. **Discard** low-quality/too-short reads (prior to all processing steps)
2. **Annotate** features (locations of adapter/primer/polyA) and identities (fusion/full-length/truncated) of reads     
3. **Trim** adapter/primer/polyA sequences  
4. **Filter** low-quality/too-short reads before output to assigned files
<br>  

## General usage
Use NanoPreP `standard` mode to get ***high-quality, non-chimeric, full-length, strand-oriented, adapter/primer-removed, polyA-removed*** reads
```
python NanoPreP \
  --mode standard \
  --output_full_length output.fq \
  --report report.json \
  input.fq
```
- `--mode standard` ← run NanoPreP with `standard` mode (see section "Modes")  
- `--output_full_length output.fq` ← output ***high-quality, non-chimeric, full-length, strand-oriented, adapter/primer-removed, polyA-removed*** reads to `output.fq`  
- `--report report.json` ← recording start/stop times, the paramters used, and the detail information of `input.fq` to `report.json`  
- `input.fq` ← input FASTQ  
<br>  

## NanoPreP output
#### TODO: why annotate reads? re-usable, time-saving, transparency, flexibility
```
@read_1 strand=0.91 full_length=1 fusion=0 ploc5=0 ploc3=0 poly5=-1 poly3=0
AGAGGCTGGCGGGAACGGGC......TTTCAAAGCCAGGCGGATTC
+
+,),+'$)'%671*%('&$%......((&'(*($%$&%&$-((84*
```
NanoPreP uses the following flags to annotate each read  
|flag|regex|default|explanation|
|:-|-|-|:-|
|`strand`|-?\d+\.\d*|0|0: unknown, > 0: sense, < 0: antisense|
|`full_length`|[0\|1]|0|0: non-full-length, 1: full-length|
|`fusion`|[0\|1]|0|0: non-chimeric/-fusion, 1: chimeric/fusion|
|`ploc5`|-?\d+|-1|-1: unknown, 0: removed, > 0: 5' adapter/primer location|
|`ploc3`|-?\d+|-1|-1: unknown, 0: removed, > 0: 3' adapter/primer location|
|`poly5`|-?\d+|-1|-1: unknown, 0: removed, > 0: 5' polymer length|
|`poly3`|-?\d+|-1|-1: unknown, 0: removed, > 0: 3' polymer length|
<br>

## Modes  
NanoPreP provides several presets for different usages  
**`standard`**: output ***high-quality, non-chimeric, full-length, strand-oriented, adapter/primer-removed, polyA-removed*** reads  
   ```
   python NanoPreP \
    --mode standard \
    --output_full_length output.fq \
    --report report.json \
    input.fq
   ```
**`annotate`**: annotate without discarding/trimming/filtering reads
   ```
   python NanoPreP --mode annotate input.fq > annotated.fq
   ```
**`report`**: generate report.json from NanoPreP-annotated FASTQ files
   ```
   python NanoPreP --mode report --report report.json annotated.fq 
   ```

Options used in each mode
|mode|options|
|:-|:-|
|`standard`|```--discard_lowq 7 ```<br> ```--p5_sense TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG``` <br>```--p3_sense A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC``` <br>```--isl5 0 130``` <br>```--isl3 -60 -1``` <br>```--pid_isl 0.7``` <br>```--pid_body 0.7``` <br>```--poly_w 6``` <br>```--poly_k 4``` <br>```--filter_short 1``` <br>```--trim_adapter``` <br>```--trim_poly``` <br>```--orientation 1``` <br>```--output_full_length output.fq``` <br>```--report report.json```|
|`annotate`|```--p5_sense TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG``` <br>```--p3_sense A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC``` <br>```--isl5 0 130``` <br>```--isl3 -60 -1``` <br>```--pid_isl 0.7``` <br>```--pid_body 0.7``` <br>```--poly_w 6``` <br>```--poly_k 4``` ```--orientation 0``` <br>```--output_fusion -``` <br>```--output_truncated -``` <br>```--output_full_length -``` <br>```--report report.json```|
|`report`|```--disable_annot```<br>```--report report.json```|
<br>



## Options
```
usage: NanoPreP [-h] [--discard_lowq int] [--discard_short int]
                [--disable_annot] [--p5_sense str] [--p3_sense str]
                [--isl5 int int] [--isl3 int int] [--pid_isl float]
                [--pid_body float] [--poly_w int] [--poly_k int]
                [--orientation int] [--trim_adapter] [--trim_poly]
                [--filter_short int] [--filter_lowq float]
                [--mode [strandard|untrimmed|polyA|annotate]]
                [--suffix_filtered str] [--report PATH] [--config PATH]
                [--output_fusion PATH] [--output_truncated PATH]
                [--output_full_length PATH]
                input.fq

positional arguments:
  input.fq              input FASTQ  

annotation options:
  --disable_annot       use this flag to disable annotation
  --p5_sense str        5' sense adatper/primer/polymer sequences
  --p3_sense str        3' sense adatper/primer/polymer sequences
  --isl5 int int        ideal searching location for 5' adapter/primer
                        sequences (e.g. 1 130)
  --isl3 int int        ideal searching location for 3' adapter/primer
                        sequences (e.g. -60 -1)
  --pid_isl float       adapter/primer percent identity cutoff (in ISL)
  --pid_body float      adapter/primer percent identity cutoff (on read body)
  --poly_w int          window size for homopolymer identification
  --poly_k int          number of monomers to be expected in the window

processing (orient/trim/filter) options:
  --discard_lowq int    discard low-quality reads prior to all processing
                        steps (default: -1)
  --discard_short int   discard too-short reads prior to all processing steps
                        (default: -1)
  --orientation int     re-orient reads (0: generic (default), 1: sense, -1:
                        antisense)
  --trim_adapter        use this flag to trim adatper/primer sequences
  --trim_poly           use this flag to trim homopolymers
  --filter_lowq float   filter low-quality reads after all trimming steps
                        (default: -1)
  --filter_short int    filter too short reads after all trimming steps
                        (default: -1)

general options:
  -h, --help            show this help message and exit
  --mode [strandard|polyA|untrimmed]
                        use parameter presets (can be overriden by `config`
                        and command line arguments)
  --config PATH         use the parameters in this config file (JSON)(can be
                        overriden by command line arguments)
  --report PATH         output report file (JSON)
  --output_fusion PATH  output fusion/chimeric reads to this file (use '-'
                        to output to stdout)
  --output_truncated PATH
                        output truncated/non-full-length reads to this file (use '-' to output to stdout)
  --output_full_length PATH
                        output full-length reads to this file (use '-' to output to stdout)
  --suffix_filtered str
                        output filtered reads with the file suffix
```
