#!/usr/bin/env python

import os
import sys

o = open(sys.argv[1], "w+")


o.write("<html> <body> <h1> InterProScan result summary page </h1> <ul>")

for filename in [
    f for f in os.listdir(sys.argv[2]) if os.path.isfile(os.path.join(sys.argv[2], f))
]:
    o.write(
        '<li><a href="%s"> %s </a></li>' % (filename, os.path.splitext(filename)[0])
    )

o.write("</ul></body></html>")
o.close()
