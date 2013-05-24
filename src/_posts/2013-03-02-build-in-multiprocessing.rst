---
title: Divide & Conquer (mutiprocessing)
layout: post
---

Galaxy has some initial support for build-in multiprocessing.
In the datatype definition you can define merge() and split() functions 
to devide & conquer your input datasets. After processing every splitted file all results 
will be merged automatically.


I commited today an extended definition of SMILES, InChI and SDF datatypes and
enabled the multiprocessing feature in a few tools, like QED and Converters.
Depending on your galaxy configuration and the available computer cores your calculations
can be X times faster now.

During processing 50.000.000 SMILES I found a small bug in the concatination routine of galaxy.
Using cat to merge 1000 of files seems not to be suitable, because the shell has a commandline limit
defined in ARG_MAX:


.. code-block:: console

   $ getconf ARG_MAX
   $ 2097152


The solution is to use python's shutil module and iterate over all splitted files.


.. code-block:: python
   :caption: Python file merging:

   for fsrc in split_files:
      shutil.copyfileobj( open( fsrc, 'rb' ), fdst )


That solution is not slower than cat, pythonic and it works for a unlimited amount of files.
Patch and Pull request is submitted as `#141`_.


.. _#141: https://bitbucket.org/galaxy/galaxy-central/pull-request/141/


