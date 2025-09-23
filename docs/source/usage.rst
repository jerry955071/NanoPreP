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

1. :code:`--beta 0.2`: :math:`F_{\beta=0.2}` optimization for adapter/primer alignment cutoffs (see :ref:`here <f_beta_optimization>`)
2. :code:`--output_full_length output.fq`: identify full-length reads (see :ref:`here <read_classification>`)
3. :code:`--trim_adapter`: trim adapter/primer (see :ref:`here <trim_ap>`)
4. :code:`--trim_poly`: trim poly A/T (see :ref:`here <trim_poly>`)
5. :code:`--orientation 1`: reorientation (see :ref:`here <reorient>`)
6. :code:`--filter_lowq 7`: filter low-quality reads (see :ref:`here <read_filter>`)



Pre-processing pipeline
-------------

.. _f_beta_optimization:

:math:`F_{\beta}` optimization
~~~~~~~~~~~~~
Identifying adapter/primer is the basis of pre-processing by NanoPrePro.
NanoPrePro optimizes adapter/primer alignment cutoff by:

1. Simulating true and random alignment
2. Find the cutoff values that best separates true and random alignments


.. _read_classification:

Identify full-length/truncated/chimeric reads
~~~~~~~~~~~~~

.. _trim_ap:

Trim adapter/primer
~~~~~~~~~~~~~

.. _trim_poly:

Trim poly A/T
~~~~~~~~~~~~~

.. _reorient:

Reorientation
~~~~~~~~~~~~~

.. _read_filter:

Filter low-quality reads
~~~~~~~~~~~~~

