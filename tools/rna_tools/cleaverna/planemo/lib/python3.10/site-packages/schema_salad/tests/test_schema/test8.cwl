#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
inputs:
  foo: string
outputs:
  bar: string
steps:
  step1:
    scatterMethod: abc
    in: []
    out: [out]
