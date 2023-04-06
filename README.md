# NanoPreP: a fully-equipped, fast, and memory-efficient pre-processor for ONT transcriptomic data

## Requirements
* Python (>= 3.7)  
* edlib (>=1.3.8) for Python

## Getting started
**Option 1:** use git clone and run NanoPreP as a python module without installation
```
git clone https://github.com/Woodformation1136/NanoPreP.git
cd NanoPreP
python -m NanoPreP --help
```
**Option 2:** use pip install
```
pip install nanoprep-ffm
nanoprep --help
```

## How NanoPreP works
### **Read annotation**  
NanoPreP first annotates **locations of adapter/primer and polyA/T** on each input read using a **"local-searching"** method, and classifies each read as either **fusion, full-length or truncated** based on the adapter/primer locations.  

For the **"local-searching"** method, adapter/primer sequences (see section [How to specify adapter/primer and polyA/T sequences](#HOWTO)) are aligned to both (1) the **"ideal searching locations"** for adapter/primer (`isl5` and `isl3`) and (2) the **"read body"** between 5' and 3' ideal searching locations.

Based on the alignment hits of adapter/primer, each read can be classified into one of the three categories:  
1. **Fusion/chimeric**: reads with a adapter/primer hit on the **"read body"**
2. **Full-length**: reads that are not **fusion/chimeric** and possess both **5' and 3' adapters/primers** on **5' and 3' "ideal searching locations"**, respecitvely.
3. **Truncated**: reads that are not **fusion/chimeric** and not **full-length**.  

During the annotation step, user can decide to skip low-quality or too-short reads using options `--skip_lowq` and `--skip_short` with the desired cutoff values, respectively. The skipped reads will not be annotated and will not be included in the final output file.

### **Read processing: trimming, filtering and orientation**
After the annotation steps, **trimming**, **filtering** and **orientation** can be performed on each read. 
- **Trimming** of adapter/primer and polyA/T sequences can be applied with the flags `--trim_adapter` and `--trim_poly`.  
- **Filtering** of low-quality or too-short sequences (after trimming) can be performed using options `--filter_lowq` and `--filter_short` with desired cutoff values.   
- **Orientation** of read strand to sense or antisense can be performed using option `--orientation`.

### **Read output**
For the output options, user can choose to output **fusion/chimeric**, **full-length**, and/or **truncated** reads using options `--output_fusion`, `--output_full_length`, and `--output_truncated`, respectively, with assigned file names (or '-' to write to stdout).

## General usage
Use NanoPreP with the `standard` mode to get ***high-quality, non-chimeric, full-length, strand-oriented, adapter/primer-removed, polyA-removed*** reads:
```
nanoprep \
  --mode standard \
  --output_full_length output.fq \
  --report report.json \
  input.fq
```
- `--mode standard` ← run NanoPreP with `standard` mode (see section [Modes](#Modes))  
- `--output_full_length output.fq` ← output full-length reads to `output.fq`  
- `--report report.json` ← output details of the run to `report.json`  
- `input.fq` ← input FASTQ  

<!-- TODO: why annotate reads? re-usable, time-saving, transparency, flexibility -->
After running this command, two output files `output.fq` and `report.json` will appear in your working directory.  The `report.json` records start/stop times, the parameters used, and the detail information of the input FASTQ file. The `output.fq` contains full-length reads processed by NanoPreP. For each processed read, NanoPreP appends the annotations to the ID line (the line started with @): 
```
@read_1 strand=0.91 full_length=1 fusion=0 ploc5=0 ploc3=0 poly5=-1 poly3=0
AGAGGCTGGCGGGAACGGGC......TTTCAAAGCCAGGCGGATTC
+
+,),+'$)'%671*%('&$%......((&'(*($%$&%&$-((84*
```
As shown in the example above, several flags are used for the annotation: 
|flag|regex|default|explanation|
|:-|-|-|:-|
|`strand`|-?\d+\.\d*|0|0: unknown; > 0: sense; < 0: antisense|
|`full_length`|[0\|1]|0|0: non-full-length; 1: full-length|
|`fusion`|[0\|1]|0|0: non-chimeric/-fusion; 1: chimeric/fusion|
|`ploc5`|-?\d+|-1|-1: unknown; 0: removed; > 0: 5' adapter/primer location|
|`ploc3`|-?\d+|-1|-1: unknown; 0: removed; > 0: 3' adapter/primer location|
|`poly5`|-?\d+|-1|-1: unknown; 0: removed; > 0: 5' polymer length|
|`poly3`|-?\d+|-1|-1: unknown; 0: removed; > 0: 3' polymer length|

According to the annotation, the example "read1" is a sense strand (`strand=0.91`), full-length (`full_length=1`), non-chimeric (`fusion=0`),  adapter/primer removed (`ploc5=0 ploc3=0`), and polyA removed (`poly3=0`) read.

## Modes <a id="Modes"></a>
In addition to the `standard` mode, NanoPreP also provides other `mode` options for different usages:  
1. **`standard`**: output ***high-quality, non-chimeric, full-length, strand-oriented, adapter/primer-removed, polyA-removed*** reads  
   ```
   nanoprep \
    --mode standard \
    --output_full_length output.fq \
    --report report.json \
    input.fq
   ```
2. **`annotate`**: annotate without skipping/trimming/filtering/orienting reads
   ```
   nanoprep --mode annotate input.fq > annotated.fq
   ```
3. **`report`**: generate report.json from NanoPreP-annotated FASTQ files
   ```
   nanoprep --mode report --report report.json annotated.fq 
   ```

Each `mode` option applies multiple options at the same time, which can be detailed as follows:
|mode|options|
|:-|:-|
|`standard`|```--discard_lowq 7 ```<br> ```--p5_sense TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG``` <br>```--p3_sense A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC``` <br>```--isl5 0 130``` <br>```--isl3 -60 -1``` <br>```--pid_isl 0.7``` <br>```--pid_body 0.7``` <br>```--poly_w 6``` <br>```--poly_k 4``` <br>```--filter_short 1``` <br>```--trim_adapter``` <br>```--trim_poly``` <br>```--orientation 1``` <br>```--output_full_length output.fq``` <br>```--report report.json```|
|`annotate`|```--p5_sense TCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTGGG``` <br>```--p3_sense A{100}GAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAAC``` <br>```--isl5 0 130``` <br>```--isl3 -60 -1``` <br>```--pid_isl 0.7``` <br>```--pid_body 0.7``` <br>```--poly_w 6``` <br>```--poly_k 4``` <br>```--orientation 0``` <br>```--output_fusion -``` <br>```--output_truncated -``` <br>```--output_full_length -``` <br>```--report report.json```|
|`report`|```--disable_annot```<br>```--report report.json```|
<br>

## How to specify adapter/primer and polyA/T sequences <a id="HOWTO"></a>
To perform read annotation, **adapter/primer** and **polyA/T** (if provided) sequences on the **sense strand** can be specified using options `--p5_sense` and `--p3_sense`. 

For example, the following command means that the 5' and 3' adatper/primer sequences on the sense strand are 'CATTC' and 'GACTA', respectively.
```
--p5_sense CATTC --p3_sense GACTA
```
If users wish to detect polyA/T tails, a pattern `N{M}` can be used to specify the location and length of polyA/T tails:
```
--p5_sense CATTC --p3_sense A{50}GACTA
```
The command above means that there are poly`"A"` tails of maximum length `"50"` bp next to the 3' adapters/primers.


## Full usage
```
usage: NanoPreP [-h] [--discard_lowq int] [--discard_short int]
                [--disable_annot] [--p5_sense str] [--p3_sense str]
                [--isl5 int int] [--isl3 int int] [--pid_isl float]
                [--pid_body float] [--poly_w int] [--poly_k int]
                [--orientation int] [--trim_adapter] [--trim_poly]
                [--filter_short int] [--filter_lowq float]
                [--mode [strandard|annotate|report]]
                [--suffix_filtered str] [--report PATH] [--config PATH]
                [--output_fusion PATH] [--output_truncated PATH]
                [--output_full_length PATH]
                input.fq

positional arguments:
  input.fq              input FASTQ  

general options:
  -h, --help            show this help message and exit
  --mode [strandard|annotate|report]
                        use parameter presets (can be overriden by `config`
                        and command line arguments)
  --config PATH         use the parameters in this config file (JSON)(can be
                        overriden by command line arguments)
  --report PATH         output report file (JSON)

annotation options:
  --skip_lowq float     skip the annotation of low-quality reads (default: -1)
  --skip_short int      skip the annotation of too-short reads (default: -1)
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

processing (trimming/filtering/orientation) options:
  --trim_adapter        use this flag to trim adatper/primer sequences
  --trim_poly           use this flag to trim homopolymers
  --filter_lowq float   filter low-quality reads after all trimming steps
                        (default: -1)
  --filter_short int    filter too short reads after all trimming steps
                        (default: -1)
  --orientation int     re-orient reads (0: generic (default), 1: sense, -1:
                        antisense)
                        
output options:
  --output_fusion PATH  output fusion/chimeric reads to this file (use '-' for
                        stdout)
  --output_truncated PATH
                        output truncated/non-full-length reads to this file
                        (use '-' for stdout)
  --output_full_length PATH
                        output full-length reads to this file (use '-' for
                        stdout)
  --suffix_filtered str
                        output filtered reads with the suffix 
```
