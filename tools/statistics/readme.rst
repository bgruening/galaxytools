Galaxy wrapper for statistical hypothesis testing with scipy
============================================================

Computes a large number of probability distributions as well as a statistical functions of any kind.
For more informations have a look at the `SciPy site`_.

.. _`SciPy site`: http://docs.scipy.org/doc/scipy/reference/stats.html


Features
========

- **Two analysis scopes**: compute tests *per row* (each row is a set of
  observations, results are appended to the row) or *per column* (each
  selected column is a set of observations, one result per column).
- **Input flexibility**: optional header row, flexible NaN handling
  (``propagate``, ``omit``, ``raise``), and per-test column selection.
- **Rich output**: configurable labels, results-only mode, and optional
  p-value correction (Bonferroni, Benjamini-Hochberg, Benjamini-Yekutieli,
  Holm, Hochberg, FDR) across the tests of a row/column.
- **Modern test coverage**: includes bootstrap and permutation tests
  (with confidence intervals), t-test confidence intervals,
  Shapiro-Wilk, Bartlett, Levene, Ansari-Bradley, Jarque-Bera,
  Cramér-von Mises, Anderson-Darling, one- and two-sample
  chisquare/power_divergence tests, and more.
- **Tool selection diagram**: see ``static/images/statistics_tool_selection.png``
  for a visual guide on picking the right test.


============
Installation
============

Should be done via the Galaxy `Tool Shed`_.
Install the following repository: https://toolshed.g2.bx.psu.edu/view/bgruening/statistical_hypothesis_testing

.. _`Tool Shed`: http://wiki.galaxyproject.org/Tool%20Shed


=======
History
=======

  - v0.1: initial release
  - v0.2: add a lot more statistics
  - v0.3: per-row/per-column scopes, header and NaN handling, output
    labels and p-value correction, bootstrap and permutation tests,
    t-test confidence intervals, updated to scipy 1.16, tool selection
    diagram




Wrapper Licence (MIT/BSD style)
===============================

Copyright (c) 2013-2015

 * Björn Gruening (bjoern dot gruening <at> gmail dot com)
 * Hui Li (lihui900116 <at> gmail dot com)

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

