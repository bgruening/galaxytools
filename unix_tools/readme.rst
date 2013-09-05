These are Galaxy wrappers for common unix text-processing tools
===============================================================

The initial work was done by Assaf Gordon and Greg Hannon's lab ( http://hannonlab.cshl.edu ) 
in Cold Spring Harbor Laboratory ( http://www.cshl.edu ).


The tools are:

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

Few improvements over the standard tools:

  * EasyJoin - A Join tool that does not require pre-sorted the files ( https://github.com/agordon/filo/blob/scripts/src/scripts/easyjoin )
  * Multi-Join - Join multiple (>2) files ( https://github.com/agordon/filo/blob/scripts/src/scripts/multijoin )
  * Find_and_Replace - Find/Replace text in a line or specific column.
  * Grep with Perl syntax - uses grep with Perl-Compatible regular expressions.
  * HTML'd Grep - grep text in a file, and produced high-lighted HTML output, for easier viewing ( uses https://github.com/agordon/filo/blob/scripts/src/scripts/sort-header )


Requirements
------------

1. Coreutils vesion 8.19 or later.
2. AWK version 4.0.1 or later.
3. SED version 4.2 *with* a special patch
4. Grep with PCRE support

These will be installed automatically with the Galaxy Tool Shed.


-------------------
NOTE About Security
-------------------

The included tools are secure (barring unintentional bugs):
The main concern might be executing system commands with awk's "system" and sed's "e" commands,
or reading/writing arbitrary files with awk's redirection and sed's "r/w" commands.
These commands are DISABLED using the "--sandbox" parameter to awk and sed.

User trying to run an awk program similar to:
 BEGIN { system("ls") }
Will get an error (in Galaxy) saying:
 fatal: 'system' function not allowed in sandbox mode.

User trying to run a SED program similar to:
 1els
will get an error (in Galaxy) saying:
 sed: -e expression #1, char 2: e/r/w commands disabled in sandbox mode

That being said, if you do find some vulnerability in these tools, please let me know and I'll try fix them.

------------
Installation
------------

Should be done with the Galaxy `Tool Shed`_.

.. _`Tool Shed`: http://wiki.galaxyproject.org/Tool%20Shed


----
TODO
----

- unit-tests
- uniqu will get a new --group funciton with the 8.22 release, its currently commended out
- also shuf will get a major improved performance with large files http://git.savannah.gnu.org/gitweb/?p=coreutils.git;a=commit;h=20d7bce0f7e57d9a98f0ee811e31c757e9fedfff
  we can remove the random feature from sort and use shuf instead
- move some advanced settings under a conditional, for example the cut tools offers to cut bytes
- cut wrapper has some output conditional magic for interval files, that needs to be checked
- comm wrapper, see the Galaxy default one
- evaluate the join wrappers against the Galaxy ones, maybe we should drop them





