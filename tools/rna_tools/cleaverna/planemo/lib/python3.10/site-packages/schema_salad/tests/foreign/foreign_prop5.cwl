#!/usr/bin/env cwl-runner
doc: |
  Foreign properties test.  The property is supposed to cross
  reference another concept node, but that node doesn't exist.
  If strict_foreign_properties is true, it should fail.  If
  false, should pass with a warning.
cwlVersion: v1.0
$schemas:
  - ../EDAM.owl
$namespaces:
  edam: http://edamontology.org/
class: CommandLineTool
inputs: []
outputs: []
baseCommand: echo
edam:has_topic: abc
