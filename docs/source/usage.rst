Usage
=====

Quick Start
-----------

Run a standard preprocessing pipeline using NanoPrePro as follows:

.. code-block:: bash

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

This command performs the following preprocessing steps and generates a report file (:code:`report.html`):

1. :code:`--beta 0.2`: performs :math:`F_{\beta=0.2}` optimization for adapter/primer alignment cutoffs (see :ref:`here <f_beta_optimization>`).
2. :code:`--output_full_length output.fq`: identifies full-length reads (see :ref:`here <read_classification>`).
3. :code:`--trim_adapter`: trims adapter/primer sequences (see :ref:`here <trim_ap>`).
4. :code:`--trim_poly`: trims poly(A/T) sequences (see :ref:`here <trim_poly>`).
5. :code:`--orientation 1`: reorients reads (see :ref:`here <reorient>`).
6. :code:`--filter_lowq 7`: filters low-quality reads (see :ref:`here <read_filter>`).

Preprocessing Pipeline
----------------------

.. _f_beta_optimization:

Step 1. :math:`F_{\beta}` Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NanoPrePro optimizes adapter/primer alignment cutoffs by:

1. Simulating both true and random alignments.  
2. Identifying cutoff values that best separate true from random alignments.  

Adapter/primer sequences provided by the user are aligned twice to each read. 
(:code:`--p5_sense 5_PRIMER_SEQUENCE` and :code:`--p3_sense A{100}3_PRIMER_SEQUENCE`)

.. note::

   :code:`A{100}` indicates that up to 100 consecutive :code:`A` nucleotides 
   may occur adjacent to the 3′ adapter/primer. These bases are **NOT** used 
   for alignment. See :ref:`Poly A/T trimming <trim_poly>` for details.

NanoPrePro then determines the alignment cutoffs that maximize the :math:`F_{\beta}` score (:code:`--beta <float>`),  
the weighted harmonic mean of precision and recall:

.. math::
   :nowrap:

   \text{Let:} \\
   \mathrm{TP} = \text{True Positives (true alignments passes the cutoffs)} \\
   \mathrm{FP} = \text{False Positives (random alignments passes the cutoffs)} \\
   \mathrm{FN} = \text{False Negatives (true alignments rejected by the cutoffs)}

.. math::
   :nowrap:

   \mathrm{precision} = \frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FP}}

.. math::
   :nowrap:
   
   \mathrm{recall} = \frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FN}}

The :math:`\beta` parameter controls the weighting of precision versus recall:

- Higher :math:`\beta` values emphasize recall.  
- Lower :math:`\beta` values emphasize precision.  

For recommended :math:`\beta` ranges for ONT datasets with different kits and chemistries,  
please refer to our :ref:`manuscript <#TODO>`.

.. math::

   F_{\beta} = (1 + \beta^2) \cdot \frac{\mathrm{precision} \cdot \mathrm{recall}}
   {(\beta^2 \cdot \mathrm{precision}) + \mathrm{recall}}

The cutoff values achieving the highest :math:`F_{\beta}` score are used for adapter/primer identification.

.. _read_classification:

Step 2. Full-Length / Truncated / Chimeric Read Classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reads are classified into three categories based on adapter/primer alignment results:

- **Full-length**: 5' and 3' adapter/primer present, no internal adapters/primers.  
- **Chimeric**: contains internal adapter/primer sequences.  
- **Truncated**: neither chimeric nor full-length.

Output files for each read type can be specified as:

- Full-length: :code:`--output_full_length` (default to standard output).  
- Chimeric: :code:`--output_fusion`.  
- Truncated: :code:`--output_truncated`.

.. _trim_ap:

Step 3. Adapter/Primer Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is activated with :code:`--trim_adapter`.  
It removes flanking (5' and/or 3') adapter/primer sequences from the output reads.

.. note::

   Trimming is applied to all requested output reads, regardless of read type.

.. _trim_poly:

Step 4. Poly(A/T) Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is activated with :code:`--trim_poly`.  
The expected length, location, and nucleotide of mono-polymers are assigned along with the primer sequence.

Use a pattern like :code:`N{M}` to specify the location and length of polyA/T tails. For example, this command tells NanoPrePro that poly :code:`A` tails of up to :code:`50` nucleotides occur adjacent to the 3' adapters/primers:

.. code::

   --p3_sense A{50}GACTA

.. note::

   Poly(A/T) trimming is applicable only if adapters/primers are trimmed. 
   Similar to adapter/primer trimming, this step can be performed on all classes of output reads. 

.. _reorient:

Step 5. Read Reorientation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Read strands are determined based on the orientation of aligned adapter/primer sequences.  
Adapter/primer sequences should be provided in the sense direction.  
Reads aligned to the reverse complement are classified as antisense.

Reorientation can be performed using :code:`--orientation [0, 1, -1]`:

- `1`: sense direction  
- `-1`: antisense  
- `0`: do not reorient

.. _read_filter:

Step 6. Filtering Low-Quality Reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Average Q-scores are calculated after trimming adapter/primer/polyA(T) sequences (if applied).  
Trimming removes low-quality regions at read termini, providing a more accurate measure of read quality.

Step 7. Output
~~~~~~~~~~~~~~

NanoPrePro produces:

- **FASTQ**: processed reads  
- **HTML report**: summary of preprocessing statistics

**FASTQ Files**  

Processed reads are saved separately for full-length, truncated, and chimeric reads.  
Output file names can be assigned with :code:`--output_full_length`, :code:`--output_truncated`, and :code:`--output_fusion`.

.. note::

   Gzip-compressed FASTQ files are supported. For example:  
   :code:`--output_full_length output.fq.gz`

Per-read annotations are appended to FASTQ read IDs.  
See :ref:`output<per_read_annotation>` for details.

**HTML Report**  

Written to the file specified by :code:`--report`.  
The report includes Q-score distributions, the proportion of full-length/truncated/chimeric reads, and adapter/primer alignment results from :math:`F_{\beta}` optimization.

The simulated alignment results help users manually picking cutoffs. 
See :ref:`output<guideline>` for guidelines on manually selecting alignment cutoffs based on simulated alignment data.
