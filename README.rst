=====
moose
=====

A tool for grouping and summarizing large alignment or kmer based
classifications into something more concise and readable.

Moose stands for nothing.  It is the name of my cat.

.. contents:: Table of Contents

authors
=======

* `Chris Rosenthal <crosenth@gmail.com>`_

about
=====

Moose groups pairwise alignments or kmer counts by taxonomy and alignment 
scores.  It works safely with large data sets utlizing the Python Data 
Analysis Library.

dependencies
============

* Python 3.x
* `Pandas <https://pandas.pydata.org/>`_ >= 0.24.0

installation
==============

moose can be installed in a few ways::

  % pip3 install moose

For developers::

  % pip3 install git://github.com/crosenth/moose.git
  # or
  % git clone git://github.com/crosenth/moose.git 
  % cd medirect
  % python3 setup.py install
