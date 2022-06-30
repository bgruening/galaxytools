#!/usr/bin/env python
import json
import argparse


output = {
    "tool": {
        "name": "hAMRonize",
        "version": "?",
        "description": "Tool for combining results",
        "configuration": {"verbose": "true", "multisetting": ["first", "second"]},
    },
    "records": [],
}

with open("input.json", "r") as handle:
    value = handle.read()

    x = json.loads(value)

    for i in x:

        start = i["input_gene_start"]
        end = i["input_gene_stop"]
        label = i["gene_name"]
        seqid = i["input_sequence_id"]

        score = i["sequence_identity"]
        if start <= end:

            record = {
                "name": seqid,
                "subregions": [
                    {
                        "start": (start - 1),
                        "end": end,
                        "label": label,
                        "details": {"score": score},
                    }
                ],
            }
            output["records"].append(record)
        else:
            record = {
                "name": seqid,
                "subregions": [
                    {
                        "start": (end - 1),
                        "end": start,
                        "label": label,
                        "details": {"score": score},
                    }
                ],
            }
            output["records"].append(record)


with open("output.json", "w") as handle:

    a = json.dumps(output, indent=2)
    handle.write(a)
