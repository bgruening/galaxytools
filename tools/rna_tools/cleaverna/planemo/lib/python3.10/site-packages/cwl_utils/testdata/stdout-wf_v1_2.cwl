#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
inputs:
  letters0:
    type: string
    default: "a0"

outputs:
  all:
    type: File
    outputSource: echo_w/txt

steps:
  echo_w:
    run: echo_v1_2.cwl
    in:
      echo_in: letters0
    out: [txt]

