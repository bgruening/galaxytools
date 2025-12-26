#!/usr/bin/env cwl-runner
doc: |
  Foreign properties test.  This property is not part of the supplied
  ontology.  If strict_foreign_properties is true, it should fail.  If
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
'edam:fake_property': abc
