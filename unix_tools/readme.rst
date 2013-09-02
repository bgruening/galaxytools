These are Galaxy wrappers for common unix text-processing tools.

Source:
http://hannonlab.cshl.edu/galaxy_unix_tools/index.html

Contact: gordon at cshl dot edu

NOTE: You must install some programs manually. See below for details.

The tools are:

* awk - The AWK programmning language ( http://www.gnu.org/software/gawk/ )
* sed - Stream Editor ( http://sed.sf.net )
* grep - Search files ( http://www.gnu.org/software/grep/ )
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
  * Sort-Header - Sort a file, while maintaining the first line as header line ( https://github.com/agordon/filo/blob/scripts/src/scripts/sort-header )
  * Find_and_Replace - Find/Replace text in a line or specific column.
  * Grep with Perl syntax - uses grep with Perl-Compatible regular expressions.
  * HTML'd Grep - grep text in a file, and produced high-lighted HTML output, for easier viewing ( uses https://github.com/agordon/filo/blob/scripts/src/scripts/sort-header )


Requirements
============
1. Coreutils vesion 8.19 or later.
2. AWK version 4.0.1 or later.
3. SED version 4.2 *with* a special patch
4. Grep with PCRE support


NOTE About Security
===================
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

That being said, if you do find some vulnerability in these tools, please let me know and I'll fix them.


Installation
============

Should be done with the Galaxy `Tool Shed`_.

.. _`Tool Shed`: http://wiki.galaxyproject.org/Tool%20Shed


