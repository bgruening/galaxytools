Galaxy wrapper for sed
======================

This wrapper is copyright 2012 by Björn Grüning.

This is a wrapper for the sed command line tool.

sed is a stream editor included in every unix-derived operating system.
That wrapper only uses a small subset of the sed functionality. Its only a wrapper for:

sed -r '$pattern' $input


WARNING:
========

No syntax check and sanitising will happen in that wrapper. This wrapper may harm your computer ;)
Nevertheless, i think it can be useable for some installations.


Installation
============

sed should be available on every unix derived operating system and belongs to the classical unix programms.
No further installation is requiered.

For more information have a look at the following pages:

http://en.wikipedia.org/wiki/Sed
http://www.gnu.org/software/sed/manual/sed.html

To install the wrapper copy the sed folder in the galaxy tools
folder and modify the tools_conf.xml file to make the tool available to Galaxy.
For example:

<toolbox>
    <tool file="text_manipulation/sed.xml" />
</toolbox>


History
=======

v0.1 - Initial public release


Wrapper Licence (MIT/BSD style)
===============================

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

