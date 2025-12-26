#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Example trivial wrapper for Java 9 compiler
requirements:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
hints:
  DockerRequirement:
    dockerPull: openjdk:9.0.1-11-slim
inputs:
  src:
    type: File
    inputBinding:
      position: 1
baseCommand: javac
arguments: ["-d", $(runtime.outdir)]
outputs:
  classfile:
    type: File
    outputBinding:
      glob: "*.class"

