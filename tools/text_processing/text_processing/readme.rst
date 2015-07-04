Galaxy wrappers for common unix text-processing tools
=====================================================

The initial work was done by Assaf Gordon and Greg Hannon's lab ( http://hannonlab.cshl.edu )
in Cold Spring Harbor Laboratory ( http://www.cshl.edu ). In late 2013 maintainence and
further development was taken over by Bjoern Gruening. Feel free to contribute any general purpose
text manipulation tool to this repository.


Tools:
------

  * awk - The AWK programmning language ( http://www.gnu.org/software/gawk/ )
  * sed - Stream Editor ( http://sed.sf.net )
  * grep - Search files ( http://www.gnu.org/software/grep/ )
  * sort_columns - Sorting every line according to there columns
  * GNU Coreutils programs ( http://www.gnu.org/software/coreutils/ ):

  * sort - sort files
  * join - join two files, based on common key field.
  * cut  - keep/discard fields from a file
  * unsorted_uniq - keep unique/duplicated lines in a file
  * sorted_uniq - keep unique/duplicated lines in a file
  * head - keep the first X lines in a file.
  * tail - keep the last X lines in a file.
  * unfold_column - unfold a column with multiple entities into multiple lines


Few improvements over the standard tools:
-----------------------------------------

  * EasyJoin - A Join tool that does not require pre-sorted the files ( https://github.com/agordon/filo/blob/scripts/src/scripts/easyjoin )
  * Multi-Join - Join multiple (>2) files ( https://github.com/agordon/filo/blob/scripts/src/scripts/multijoin )
  * Find_and_Replace - Find/Replace text in a line or specific column.
  * Grep with Perl syntax - uses grep with Perl-Compatible regular expressions.
  * HTML'd Grep - grep text in a file, and produced high-lighted HTML output, for easier viewing ( uses https://github.com/agordon/filo/blob/scripts/src/scripts/sort-header )


Requirements:
-------------

    * Coreutils vesion 8.22 or later.
    * AWK version 4.0.1 or later.
    * SED version 4.2 *with* a special patch
    * Grep with PCRE support

All dependencies will be installed automatically with the Galaxy `Tool Shed`_ and the following repository: https://toolshed.g2.bx.psu.edu/view/bgruening/text_processing


-------------------
NOTE About Security
-------------------

The included tools are secure (barring unintentional bugs):
The main concern might be executing system commands with awk's "system" and sed's "e" commands,
or reading/writing arbitrary files with awk's redirection and sed's "r/w" commands.
These commands are DISABLED using the "--sandbox" parameter to awk and sed.

User trying to run an awk program similar to::

  BEGIN { system("ls") }

Will get an error (in Galaxy) saying::

  fatal: 'system' function not allowed in sandbox mode.

User trying to run a SED program similar to::

  1els

will get an error (in Galaxy) saying::

  sed: -e expression #1, char 2: e/r/w commands disabled in sandbox mode

That being said, if you do find some vulnerability in these tools, please let me know and I'll try fix them.

------------
Installation
------------

Should be done via the Galaxy `Tool Shed`_.
Install the following repository: https://toolshed.g2.bx.psu.edu/view/bgruening/text_processing

.. _`Tool Shed`: http://wiki.galaxyproject.org/Tool%20Shed


----
TODO
----

 * add shuf, we can remove the random feature from sort and use shuf instead
 * move some advanced settings under a conditional, for example the cut tools offers to cut bytes
 * cut wrapper has some output conditional magic for interval files, that needs to be checked
 * comm wrapper, see the Galaxy default one
 * evaluate the join wrappers against the Galaxy ones, maybe we should drop them


-------
License
-------

  * Copyright (c) 2009-2013   A. Gordon  (gordon <at> cshl dot edu)
  * Copyright (c) 2013-2015   B. Gruening  (bjoern dot gruening <at> gmail dot com)


Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
