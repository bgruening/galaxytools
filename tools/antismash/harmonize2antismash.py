#!/usr/bin/env python
import json
import sys


output = {
    "tool": {
        "name": "hAMRonize",
        "version": "?",
        "description": "Tool for combining results",
        "configuration": {
            "verbose": "true", "multisetting": ["first", "second"]
        },
    },
    "records": [],
}

with open(sys.argv[1], "r") as handle:
    records = json.load(handle)

    for i in records:
        start = i["input_gene_start"]
        end = i["input_gene_stop"]
        label = i["gene_name"]
        seqid = i["input_sequence_id"]

        score = i["sequence_identity"]
        if start <= end:
            fstart = start - 1
            fend = end
        else:
            fstart = end - 1
            fend = start

        record = {
            "name": seqid,
            "subregions": [
                {
                    "start": fstart,
                    "end": fend,
                    "label": label,
                    "details": {"score": score},
                }
            ],
        }
        output["records"].append(record)


with open(sys.argv[2], "w") as handle:
    json.dump(output, handle, indent=2)
