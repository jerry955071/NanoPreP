Usage
=====

Quick Start
-----------

Suppose your sequence library looks like this:

.. image:: images/library_construction.png
   :alt: library construction


Run a standard pre-processing pipeline using NanoPrePro as follows:

.. code-block:: bash

   nanoprepro \
      --input_fq input.fq \
      --p5_sense ATCGATCG \ # 5' adapter/primer sequence (sense strand; 5' to 3')
      --p3_sense A{20}GCAATGA \ # 3' adapter/primer sequence (sense strand; 5' to 3')
      --beta 0.2 \
      --output_full_length output.fq \
      --trim_adapter \
      --trim_poly \
      --orientation 1 \
      --filter_lowq 7 \
      --report report.html

This command performs the following pre-processing steps and 
generates a report file (:code:`report.html`):

1. :code:`--beta 0.2`: performs :math:`F_{\beta=0.2}` optimization for adapter/primer alignment cutoffs (see :ref:`Step 1 <f_beta_optimization>`).
2. :code:`--output_full_length output.fq`: identifies full-length reads (see :ref:`Step 2 <read_classification>`).
3. :code:`--trim_adapter`: trims adapter/primer sequences (see :ref:`Step 3 <trim_ap>`).
4. :code:`--trim_poly`: trims poly(A/T) sequences (see :ref:`Step 4 <trim_poly>`).
5. :code:`--orientation 1`: reorients reads to sense strand (see :ref:`Step 5 <reorient>`).
6. :code:`--filter_lowq 7`: filters low-quality (avg. Q-score < 7) reads (see :ref:`Step 6 <read_filter>`).

pre-processing Pipeline
----------------------

.. _f_beta_optimization:

Step 1. :math:`F_{\beta}` Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NanoPrePro optimizes adapter/primer alignment cutoffs by:

1. Simulating both true and random alignments.  
2. Identifying cutoff values that best separate true from random alignments.  

First, the adapter/primer sequences provided by the user are aligned twice to each read. 
(:code:`--p5_sense ATCGATCG` and :code:`--p3_sense A{20}GCAATGA`)

NanoPrePro then search for the alignment cutoffs that maximize the :math:`F_{\beta}` score 
(:code:`--beta <float>`), the weighted harmonic mean of precision and recall:

:math:`\text{Precision} = \frac{\text{true alignments that pass the cutoffs}}{\text{true alignments that pass the cutoffs} + \text{random alignments that pass the cutoffs}}`

:math:`\text{Recall} = \frac{\text{true alignments that pass the cutoffs}}{\text{all true alignments}}`

:math:`F_{\beta}` score is calculated as:

.. math::

   F_{\beta} = (1 + \beta^2) \cdot \frac{\mathrm{precision} \cdot \mathrm{recall}}
   {(\beta^2 \cdot \mathrm{precision}) + \mathrm{recall}}

The :math:`\beta` parameter controls the weighting of precision versus recall:

- Higher :math:`\beta` values emphasize recall.  
- Lower :math:`\beta` values emphasize precision.  

The alignment cutoff values achieving the highest :math:`F_{\beta}` score are used for adapter/primer identification.

.. note::

   :code:`A{20}` indicates that up to 20 consecutive :code:`A` nucleotides 
   may occur adjacent to the 3′ adapter/primer. These bases are **NOT** used 
   for alignment. See :ref:`Poly A/T trimming <trim_poly>` for details.


.. _read_classification:

Step 2. Full-Length / Truncated / Chimeric Read Classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reads are classified into three categories based on adapter/primer alignment results:

- **Full-length**: 5' and 3' adapter/primer present, no internal adapters/primers.  
- **Chimeric**: contains internal adapter/primer sequences.  
- **Truncated**: not chimeric and not full-length.

Output files for each read type can be specified using:

- Full-length: :code:`--output_full_length` (default to standard output).  
- Chimeric: :code:`--output_fusion`.  
- Truncated: :code:`--output_truncated`.

.. _trim_ap:

Step 3. Adapter/Primer Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is activated with :code:`--trim_adapter`.  
It trims adapter/primer sequences from the output reads.

.. note::

   Trimming is applied to all requested output reads, regardless of read type.

.. _trim_poly:

Step 4. Poly(A/T) Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is activated with :code:`--trim_poly`.  
The expected length, location, and nucleotide of mono-polymers are assigned along with the primer sequence.

Use a pattern like :code:`N{M}` to specify the location and length of poly(A/T) tails. 
For example, this command tells NanoPrePro that poly :code:`A` tails of up to :code:`20` nucleotides are adjacent to the 3' adapters/primers:

.. code::

   --p3_sense A{20}GCAATGA

NanoPrePro then use a sliding window approach to identify and trim poly(A/T) sequences.
The window size is set by :code:`--poly_w <int>` (default: 6).
The minimum number of :code:`A` or :code:`T` bases in the window is set by :code:`--poly_k <int>` (default: 4).
The length of poly(A/T) tails would be recorded in the ID line of each read (see :ref:`Output Documentation<per_read_annotation>`).

.. note::

   Poly(A/T) trimming is applicable only if adapters/primers are trimmed. 
   Similar to adapter/primer trimming, this step can be performed on all classes of output reads. 

.. _reorient:

Step 5. Read Reorientation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Read strands are determined based on the orientation of aligned adapters/primers.  
Adapter/primer sequences should be provided in the sense direction (:code:`--p5_sense` , :code:`--p3_sense`).  
Reads are determined antisense if adapters/primers are aligned in the antisense direction.

Reorientation can be performed using :code:`--orientation 1/-1/0`:

- `1`: sense direction  
- `-1`: antisense  
- `0`: do not reorient

.. _read_filter:

Step 6. Filtering Low-Quality Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Average Q-scores are calculated after trimming adapter/primer/polyA(T) sequences (if applied).  
Trimming removes low-quality regions at read termini, providing a more accurate measure of read quality.
The threshold for filtering low-quality reads can be set with :code:`--filter_lowq <int>`.

Step 7. Output
~~~~~~~~~~~~~~

NanoPrePro produces:

- **FASTQ**: processed reads  
- **HTML report**: summary of pre-processing statistics

**FASTQ Files**  

Processed reads are saved separately for full-length, truncated, and chimeric reads.  
Output file names can be assigned with :code:`--output_full_length`, :code:`--output_truncated`, and :code:`--output_fusion`.

.. note::

   Gzip-compressed FASTQ files are supported. For example:  
   :code:`--output_full_length output.fq.gz`

Per-read annotations are appended to FASTQ read IDs.  
See :ref:`Output Documentation<per_read_annotation>` for details.

**HTML Report**  

Written to the file specified by :code:`--report`.  
The report includes Q-score distributions, the proportion of full-length/truncated/chimeric reads, and adapter/primer alignment results from :math:`F_{\beta}` optimization.

The simulated alignment results help users manually picking cutoffs. 
See :ref:`Output Documentation<guideline>` for guidelines on manually selecting alignment cutoffs based on simulated alignment data.
