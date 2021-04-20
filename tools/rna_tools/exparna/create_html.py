#!/usr/bin/env python

import os
import sys
import zipfile

os.chdir(sys.argv[2])
o = open("results.html", "w+")


o.write("<html> <body> <h1> ExpaRNA Result </h1>")

for filename in os.listdir(sys.argv[2]):
    if os.path.isfile(os.path.join(sys.argv[2], filename)) and False == (
        filename.endswith("epm")
        or filename.endswith("fa")
        or filename.endswith("aln")
        or filename.endswith("html")
    ):
        o.write('<img src="%s" /><br />' % (filename))

o.write("</body></html>")
o.close()

# create zip file

zf = zipfile.ZipFile(sys.argv[1], mode="w")
for files in os.listdir(sys.argv[2]):
    zf.write(files)
zf.close()
