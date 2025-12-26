#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: geoweaver workflow

inputs:
  code_folder:
    type: File

outputs: []

steps:
  run: HelloWorld.cwl
  in:
    folder: code_folder
  out: []
