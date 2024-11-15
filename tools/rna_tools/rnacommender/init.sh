#!/usr/bin/bash

if [ ! -d "AURA_Human_data" ]; then
  wget http://www.googledrive.com/host/0B9v5_ppcfmgWNTJzVjlkc0pCMVU
  tar -xf 0B9v5_ppcfmgWNTJzVjlkc0pCMVU
  rm 0B9v5_ppcfmgWNTJzVjlkc0pCMVU
fi