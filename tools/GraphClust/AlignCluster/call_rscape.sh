#! /bin/bash

num_bps=$(grep "SS_cons" $2  | awk '{print $3}' | tr -d '.' | tr -d '\n' | wc -c)
if [ $num_bps -gt 2 ]; then
    R-scape --outdir $1 $2
else
    echo "No consensus structure, R-scape skipped";
fi