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

This performs the following pre-processing steps:

1. :ref:`:math:`F_{\beta}` optimization <f_beta_optimization>` (:code:`--beta 0.2`)
2. Identify full-length/truncated/chimeric reads (:code:`--output_full_length output.fq`)
3. Trim adapter/primer (:code:`--trim_adapter`)
4. Trim poly A/T (:code:`--trim_poly`)
5. Reorientation (:code:`--orientation 1`)
6. Filter low-quality reads (:code:`--filter_lowq 7`)



Pre-processing pipeline
-------------

.. _f_beta_optimization:

:math:`F_{\beta}` optimization
~~~~~~~~~~~~~
Identifying adapter/primer is the basis of pre-processing by NanoPrePro.
NanoPrePro optimizes adapter/primer alignment cutoff by:

1. Simulating true and random alignment
2. Find the cutoff values that best separates true and random alignments



Identify full-length/truncated/chimeric reads
~~~~~~~~~~~~~

Trim adapter/primer
~~~~~~~~~~~~~

Trim poly A/T
~~~~~~~~~~~~~

Reorientation
~~~~~~~~~~~~~

Filter low-quality reads
~~~~~~~~~~~~~

