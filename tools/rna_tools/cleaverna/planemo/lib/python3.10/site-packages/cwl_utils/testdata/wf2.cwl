#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
inputs:
  in1: types/testtypes.yml#my_boolean_array
  in2: 
    type: types/testtypes.yml#my_enum

steps:
  s1:
    run: remote-cwl/wf1.cwl
    in:
      - id: in1
        source: "#in1"
    out: [out1]
  s2:
    run: remote-cwl/tool1.cwl
    in:
      in1: s1/out1
    out: [out1]

outputs:
  out1: 
    type: string 
    outputSource: "#s2/out1"

requirements:
  SchemaDefRequirement:
    types:
      - $import: types/testtypes.yml
