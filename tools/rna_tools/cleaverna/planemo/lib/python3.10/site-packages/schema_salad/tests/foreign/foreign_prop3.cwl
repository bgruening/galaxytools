#!/usr/bin/env cwl-runner
doc: |
  Foreign properties test.  Should pass since the edam namespace is declared,
  the ontology is imported, the property is valid in the ontology, and the
  node referenced in the property is also valid in the ontology.
cwlVersion: v1.0
$schemas:
  - ../EDAM.owl
$namespaces:
  edam: http://edamontology.org/
class: CommandLineTool
inputs: []
outputs: []
baseCommand: echo
edam:has_topic: edam:topic_0003
