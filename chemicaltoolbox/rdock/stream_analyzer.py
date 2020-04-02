#!/usr/bin/env python

import sys
 

error_counter = 0

for line in sys.stdin:
    if 'RBT_DOCKING_ERROR' in line:
        error_counter += 1
        if error_counter == 10:
            print('Found exit. Terminating the program')
            exit(23)
    print(line, end = '')
