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

This performs the following pre-processing steps and created a report file :code:`--report report.html`:

1. :code:`--beta 0.2`: performs :math:`F_{\beta=0.2}` optimization for adapter/primer alignment cutoffs (see :ref:`here <f_beta_optimization>`)
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
NanoPrePro optimizes adapter/primer alignment cutoff by:

1. Simulating true and random alignment
2. Find the cutoff values that best separates true and random alignments

First, the adapter/primer sequences provided by user (:code:`--p5_sense 5_PRIMER_SEQUENCE` and :code:`--p3_sense A{100}3_PRIMER_SEQUENCE`) 
are aligned twice to each reads.

.. note::

   The expression :code:`A{100}` means that there are at most 100 :code:`A`s in front of the 3' adapter/primer and is **NOT** used for alignment. 
   See :ref:`Poly A/T trimming<trim_poly>` for details.

Second, NanoPrePro finds the alignment cutoffs that maximize the :math:`F_{\beta}` score (:code:`--beta <float>`), 
the weighted harmonic mean of precision and recall, with precision and recall calculated as: 

.. math::

   \text{precision} = \frac{true alignments accepted by the cutoffs}{random alignments accepted by the cutoffs + true alignments accepted by the cutoffs}

   \text{recall} = \frac{true alignments accepted by the cutoffs}{all true alignments}

The weight on precision and recall is controlled by the :math:`\beta` value. 
Higher :math:`\beta` favors recall and lower :math:`\beta` favors precision.
For the recommended range of :math:`\beta` values for ONT datasets of varying 
kit and chemistry, please refer to our :ref:`manuscript<#TODO>`.

.. math::

   F_{\beta} = (1 + \beta^2) \cdot \frac{\text{precision} \cdot \text{recall}}
   {(\beta^2 \cdot \text{precision}) + \text{recall}}


The cutoff values achieving the highest :math:`F_{\beta}` score will be used to identify adapters/primers.

.. _read_classification:

Full-length/truncated/chimeric reads identification
~~~~~~~~~~~~~


.. _trim_ap:

Adapter/primer trimming
~~~~~~~~~~~~~

.. _trim_poly:

Poly A/T trimming
~~~~~~~~~~~~~

.. _reorient:

Re-orientation
~~~~~~~~~~~~~

.. _read_filter:

Filter low-quality reads
~~~~~~~~~~~~~

