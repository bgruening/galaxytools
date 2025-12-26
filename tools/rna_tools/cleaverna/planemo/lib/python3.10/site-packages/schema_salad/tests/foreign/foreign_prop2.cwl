#!/usr/bin/env cwl-runner
doc: |
  Foreign properties test.  If strict_foreign_properties is true, it
  should fail.  If false, should pass no warnings.
cwlVersion: v1.0
class: CommandLineTool
$namespaces:
  edam: http://edamontology.org/
inputs: []
outputs: []
baseCommand: echo
'edam:has_topic': abc
