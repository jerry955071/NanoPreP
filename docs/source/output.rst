Output
======

NanoPrePro generates:

- Pre-processed FASTQ
- HTML report


Processed FASTQ
-------------
NanoPrePro appends summaries for pre-processing to the ID lines in FASTQ records.

Pre-processing summaries were recorded using the following flags:  

+-------------+----------------+---------+----------------------------------------------------------+
| flag        | regex          | default | explanation                                              |
+=============+================+=========+==========================================================+
| ``strand``  | -?\d+\.\d*     | 0       | 0: unknown; > 0: sense; < 0: antisense                   |
+-------------+----------------+---------+----------------------------------------------------------+
| ``full_length`` | [0\|1]     | 0       | 0: non-full-length; 1: full-length                       |
+-------------+----------------+---------+----------------------------------------------------------+
| ``fusion``  | [0\|1]         | 0       | 0: non-chimeric/-fusion; 1: chimeric/fusion              |
+-------------+----------------+---------+----------------------------------------------------------+
| ``ploc5``   | -?\d+          | -1      | -1: unknown; 0: removed; > 0: 5' adapter/primer location |
+-------------+----------------+---------+----------------------------------------------------------+
| ``ploc3``   | -?\d+          | -1      | -1: unknown; 0: removed; > 0: 3' adapter/primer location |
+-------------+----------------+---------+----------------------------------------------------------+
| ``poly5``   | -?\d+          | -1      | 0: unknown; > 0: 5' polymer length; < 0: trimmed 5' poly |
+-------------+----------------+---------+----------------------------------------------------------+
| ``poly3``   | -?\d+          | -1      | 0: unknown; > 0: 3' polymer length; < 0: trimmed 3' poly |
+-------------+----------------+---------+----------------------------------------------------------+

Example:

.. code-block:: bash
    @read_1 strand=0.91 full_length=1 fusion=0 ploc5=0 ploc3=0 poly5=0 poly3=-20
    AGAGGCTGGCGGGAACGGGC......TTTCAAAGCCAGGCGGATTC
    +
    +,),+'$)'%671*%('&$%......((&'(*($%$&%&$-((84*

According to the flags, the example "read1" is a sense strand (`strand=0.91`), full-length (`full_length=1`), non-chimeric (`fusion=0`),  adapter/primer removed (`ploc5=0 ploc3=0`), and polyA removed (`poly3=-20`) read.


HTML report
-------------
- Quality score histogram
- Proportion of filtered/passed full-length/truncated/chimeric reads
- Adapter/primer alignment cutoff optimization result

Adapter/primer alignment cutoff optimization result
~~~~~~~~~~~~~
How to pick a cutoff based on this result?
