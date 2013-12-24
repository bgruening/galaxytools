#!/usr/bin/env python

import os
import sys
import zipfile


os.chdir(sys.argv[2])

# create zip file in current directory named first argument in folder second argument

zf = zipfile.ZipFile(sys.argv[1], mode='w')
for files in os.listdir(sys.argv[2]):
        zf.write(files)
zf.close()
