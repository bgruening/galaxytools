#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

inputs:
  irec:
    type:
      name: irec
      type: record
      fields:
      - name: ifoo
        type: File
      - name: ibar
        type: File

outputs:
  orec:
    type:
      name: orec
      type: record
      fields:
      - name: ofoo
        type: File
      - name: obar
        type: File
    outputSource: step1/orec

steps:
  step1:
    run: record-output_v1_2.cwl
    in:
      irec: irec
    out: [orec]