#!/usr/bin/env bash

lib=NAM.EDTA1.8.0.MTEC02052020.TElib.fa
genome=($1)

RepeatMasker -pa 36 -q -div 40 -lib $lib -cutoff 225 -gff $genome