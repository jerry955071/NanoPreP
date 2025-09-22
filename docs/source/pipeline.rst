Pipeline Workflow
=================

The NanoPrePro pipeline includes the following steps:

1. **Adapter/primer trimming**  
   Removes adapter/primer sequences of varying length.

2. **Alignment**  
   Uses minimap2 to align reads against reference.

3. **Filtering**  
   Applies length, similarity, and location cutoffs.

4. **Statistics & Visualization**  
   Generates QC plots and summary statistics.

Workflow diagram
----------------

.. image:: _static/workflow.png
   :alt: NanoPrePro workflow
   :align: center
   :width: 80%
