Galaxy datatypes for HOMER tools
================================

These HOMER datatypes are copyright 2013 by Björn Grüning.

See the licence text below.


History
=======

======= ======================================================================
Version Changes
------- ----------------------------------------------------------------------
v0.0.1  - First release.
======= ======================================================================


Installation
============

Doing this automatically via the Galaxy Tool Shed is probably simplest.


Manual Installation
===================

Normally you would install this via the Galaxy ToolShed, which would move
the provided homer.py file into a suitable location and process the
datatypes_conf.xml entry to be combined with your local configuration.

However, if you really want to this should work for a manual install. Add
the following lines to the datatypes_conf.xml file in the Galaxy main folder::

    <datatype extension="homer_tagdir" type="galaxy.datatypes.homer:TagDirectory" mimetype="text/html" display_in_upload="false"/>

Also create the file lib/galaxy/datatypes/homer.py by moving, copying or linking
the homer.py file provided in this tar-ball.  Finally add 'import homer' near
the start of file lib/galaxy/datatypes/registry.py (after the other import
lines).


Bug Reports
===========

You can file an issue here https://github.com/bgruening/galaxytools/issues or ask
us on the Galaxy development list http://lists.bx.psu.edu/listinfo/galaxy-dev


Developers
==========

Development is happening here:

    https://github.com/bgruening/galaxytools/


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

NOTE: This is the licence for the Galaxy HOMER datatypes **only**. HOMER
and associated data files are available and licenced separately.
