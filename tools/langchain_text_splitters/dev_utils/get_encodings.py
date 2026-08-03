#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from tiktoken import get_encoding

TIKTOKEN_CACHE_DIR = "/tmp/tiktoken-cache"
os.environ["TIKTOKEN_CACHE_DIR"] = TIKTOKEN_CACHE_DIR

TIKTOKEN_ENCODINGS = [
    "o200k_harmony",
    "o200k_base",
    "cl100k_base",
    "p50k_edit",
    "p50k_base",
    "r50k_base",
    "gpt2",
]

Path(TIKTOKEN_CACHE_DIR).mkdir(parents=True, exist_ok=True)

print(f"Populating tiktoken cache directory {TIKTOKEN_CACHE_DIR}")

for encoding_name in TIKTOKEN_ENCODINGS:
    print(f"Getting {encoding_name}...")
    get_encoding(encoding_name)
