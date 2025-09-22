Usage
=====

Quick start
-------------

Standard pre-processing pipeline using NanoPrePro:

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

This performs the following steps:
1. :math:`F_{\beta}` optimization for adapter/primer alignment cutoffs
2. Classify reads in to full-length/truncated/chimeric
3. Remove adapter/primer sequences
4. Remove poly A/T adjacent to adapters/primers
5. Flip reads into sense direction
6. Filter low-quality reads after trimming

Pre-processing pipeline
-------------

:math:`F_{\beta}` optimization
~~~~~~~~~~~~~
NanoPrePro optimizes adapter/primer alignment cutoff by:
(1) simulating true and random alignment
(2) find the cutoff values that best separates true and random alignments

True and random alignments were simulated by aligning adapter/primer sequence twice to each read.
`--beta 0.2`

Read classification
~~~~~~~~~~~~~

Trim adapter/primer
~~~~~~~~~~~~~

Trim poly A/T
~~~~~~~~~~~~~

Reorientation
~~~~~~~~~~~~~

Filter low-quality/short reads
~~~~~~~~~~~~~

