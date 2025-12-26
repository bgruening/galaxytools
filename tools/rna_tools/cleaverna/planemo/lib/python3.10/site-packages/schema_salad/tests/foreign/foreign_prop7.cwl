#!/usr/bin/env cwl-runner
doc: |
  Foreign properties test.  The property
  reference an unsupported mailto: URI
cwlVersion: v1.0
$schemas:
  - ../EDAM.owl
$namespaces:
  edam: http://edamontology.org/
class: CommandLineTool
inputs: []
outputs: []
baseCommand: echo
edam:has_topic: urn:blurb
