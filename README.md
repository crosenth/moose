moose
=====

A tool for grouping and summarizing large alignment or kmer based
classifications into something more concise and readable.

Moose stands for nothing.  It is the name of my cat.

.. contents:: Table of Contents

authors
=======

* [Chris Rosenthal](<crosenth@gmail.com)

about
=====

Moose groups pairwise alignments or kmer counts by taxonomy and alignment 
scores.  It works safely with large data sets utlizing the Python Data 
Analysis Library.

dependencies
============

* Python 3.x
* [Pandas](https://pandas.pydata.org/) >= 0.24.0

installation
==============

moose can be installed in a few ways::

  `% pip3 install moose`

For developers::

  ```
  % pip3 install git://github.com/crosenth/moose.git
  # or
  % git clone git://github.com/crosenth/moose.git 
  % cd medirect
  % python3 setup.py install
  ```

examples
========

For instructions on creating a local blast nt database see the NCBI
walkthrough [here](https://www.ncbi.nlm.nih.gov/sites/books/NBK537770/)
and [here](https://www.ncbi.nlm.nih.gov/sites/books/NBK279688/).

The simplest example takes a blast "10 qaccver saccver pident staxid" and
outputs a table of species level taxonomy results::

  ```
  blastn -outfmt "10 qaccver saccver pident staxid" -query sequences.fasta -out blast.csv
  wc --lines blast.csv
  ...
  classify --columns qaccver,saccver,pident,staxid blast.csv
  |----------------+---------------+-----------------------------------------------------------+-----------+-------------+-------------+---------------+-------+----------+-----------|
  | specimen       | assignment_id | assignment                                                | best_rank | max_percent | min_percent | min_threshold | reads | clusters | pct_reads |
  |----------------+---------------+-----------------------------------------------------------+-----------+-------------+-------------+---------------+-------+----------+-----------|
  | sv-048:680-16  | 0             | Bacteria*;Escherichia coli;Staphylococcus                 | species   | 100.00      | 93.10       | 0.00          | 1     | 1        | 100.00    |
  | sv-096:680-16  | 0             | Bacteroidetes*;uncultured bacterium*/organism*            | species   | 100.00      | 92.31       | 0.00          | 1     | 1        | 100.00    |
  | sv-6669:183-1  | 0             | Homo sapiens                                              | species   | 99.67       | 99.67       | 0.00          | 1     | 1        | 100.00    |
  | sv-6685:180-17 | 0             | Homo sapiens;Pan troglodytes                              | species   | 97.07       | 86.54       | 0.00          | 1     | 1        | 100.00    |
  | sv-6686:180-9  | 0             | Homo sapiens*                                             | species   | 100.00      | 89.95       | 0.00          | 1     | 1        | 100.00    |
  | sv-6687:187-13 | 0             | Cutibacterium acnes*;Propionibacterium sp. oral taxon 193 | species   | 100.00      | 89.64       | 0.00          | 1     | 1        | 100.00    |
  | sv-6688:183-1  | 0             | Homo sapiens*;eukaryotic synthetic construct*             | species   | 100.00      | 99.56       | 0.00          | 1     | 1        | 100.00    |
  | sv-6689:187-9  | 0             | Homo sapiens                                              | species   | 99.50       | 99.50       | 0.00          | 1     | 1        | 100.00    |
  | sv-6690:183-1  | 0             | Homo sapiens*;Macaca mulatta                              | species   | 100.00      | 93.33       | 0.00          | 1     | 1        | 100.00    |
  |----------------+---------------+-----------------------------------------------------------+-----------+-------------+-------------+---------------+-------+----------+-----------|
  ```
When classify is run for the first time it will generate a lineages
to help group and collapse classification results.  To save time it is a good 
idea to save the generated lineages table to save time in subsequent runs::

  ```
  classify --lineages-out lineages.csv --columns qaccver,saccver,pident,staxid blast.csv
  classify --lineages lineages.csv --columns qaccver,saccver,pident,staxid blast.csv
  ...
  ```

The next
